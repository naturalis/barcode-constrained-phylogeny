import subprocess
import shlex
import re
from typing import Dict, List, Optional
from bactria.config import Config
from bactria.logger import get_formatted_logger


class ToolRunner:
    """
    A superclass for running command-line tools with configurable parameters and intelligent logging.
    This superclass is intended for handling the common tasks of running the command-line tools
    blastn, hmmalign, raxml-ng, sqlite3, and the megatree-* tools.

    This class provides a framework for executing command-line tools, handling their output,
    and logging the results. It uses a Config object for initialization and supports dynamic
    parameter setting. The class intelligently parses tool output to categorize log messages
    by severity.

    Attributes:
        config (Config): Configuration object containing tool and logger settings.
        logger: Configured logger instance for output handling.
        tool_name (str): Name of the command-line tool to be run.
        parameters (Dict[str, str]): Dictionary of command-line parameters for the tool.
    """

    def __init__(self, config: Config) -> None:
        """
        Initialize the ToolRunner with a configuration object.

        :param config: Configuration object containing tool and logger settings.
        :type config: Config
        """
        self.config = config
        self.logger = get_formatted_logger(self.__class__.__name__, config)
        self.tool_name: str = config.get('tool_name', '')
        self.parameters: Dict[str, str] = {}
        self._load_parameters_from_config()
        self._compile_log_level_regexes()

    def _load_parameters_from_config(self) -> None:
        """
        Load tool parameters from the configuration object.
        """
        tool_params = self.config.get('tool_parameters', {})
        for key, value in tool_params.items():
            self.set_parameter(key, value)

    def _compile_log_level_regexes(self) -> None:
        """
        Compile regex patterns for identifying log levels in tool output.
        """
        self.log_level_patterns: Dict[str, re.Pattern] = {
            'debug': re.compile(r'\b(?:debug)\b', re.IGNORECASE),
            'info': re.compile(r'\b(?:info|information)\b', re.IGNORECASE),
            'warning': re.compile(r'\b(?:warn(?:ing)?)\b', re.IGNORECASE),
            'error': re.compile(r'\b(?:error|critical|fatal)\b', re.IGNORECASE)
        }

    def set_parameter(self, key: str, value: str) -> None:
        """
        Set a command-line parameter for the tool.

        :param key: Parameter name.
        :type key: str
        :param value: Parameter value.
        :type value: str
        """
        self.logger.debug(f'Setting parameter: {key} = {value}')
        self.parameters[key] = value

    def get_parameter(self, key: str, default: Optional[str] = None) -> Optional[str]:
        """
        Get a command-line parameter value.

        :param key: Parameter name.
        :type key: str
        :param default: Default value if parameter is not set.
        :type default: Optional[str]
        :return: The parameter value or the default if not set.
        :rtype: Optional[str]
        """
        return self.parameters.get(key, default)

    def build_command(self) -> List[str]:
        """
        Build the command-line command based on the tool name and parameters.

        :return: The command as a list of strings, ready for subprocess execution.
        :rtype: List[str]
        """
        command = [self.tool_name]
        for key, value in self.parameters.items():
            if value is not None:
                if len(key) == 1:
                    command.append(f"-{key}")
                else:
                    command.append(f"--{key}")
                command.append(str(value))
        return command

    def _determine_log_level(self, line: str) -> str:
        """
        Determine the appropriate log level for a given output line.

        :param line: A line of output from the tool.
        :type line: str
        :return: The determined log level ('debug', 'info', 'warning', or 'error').
        :rtype: str
        """
        for level, pattern in self.log_level_patterns.items():
            if pattern.search(line):
                return level
        return 'info'  # Default to info if no match

    def _log_output(self, line: str, stream: str) -> None:
        """
        Log a line of output with the appropriate log level.

        :param line: A line of output from the tool.
        :type line: str
        :param stream: The stream the output came from ('stdout' or 'stderr').
        :type stream: str
        """
        line = line.strip()
        if not line:
            return

        if stream == 'stdout':
            self.logger.info(line)
        elif stream == 'stderr':
            log_level = self._determine_log_level(line)
            getattr(self.logger, log_level)(line)

    def run(self) -> int:
        """
        Run the command-line tool and handle its output.

        This method executes the tool, captures its output, determines appropriate log levels,
        and logs the output accordingly. It also handles any exceptions that occur during execution.

        :return: The return code of the command-line tool.
        :rtype: int
        :raises Exception: If an error occurs during command execution.
        """
        command = self.build_command()
        command_str = ' '.join(shlex.quote(str(arg)) for arg in command)
        self.logger.info(f"Running command: {command_str}")

        try:
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
                universal_newlines=True
            )

            # Read output until the process terminates
            for line in process.stdout:
                self._log_output(line, 'stdout')
            for line in process.stderr:
                self._log_output(line, 'stderr')

            # Wait for the process to complete and get the return code
            return_code = process.wait()

            if return_code != 0:
                self.logger.error(f"Command failed with return code {return_code}")
            else:
                self.logger.info("Command completed successfully")

            return return_code

        except Exception as e:
            self.logger.error(f"Error running command: {str(e)}")
            raise
