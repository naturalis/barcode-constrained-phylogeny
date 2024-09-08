from typing import List
from bactria.tool_runner import ToolRunner
from bactria.config import Config
import logging


class MegatreeLoader(ToolRunner):
    """
    A subclass of ToolRunner specifically for running megatree-loader.

    Examples:
        >>> config = Config()
        >>> config.load_config('path/to/config.yaml')
        >>> megatree_runner = MegatreeLoaderRunner(config)
        >>> megatree_runner.set_infile('input_tree.newick')
        >>> megatree_runner.set_dbfile('output_db.sqlite')
        >>> return_code = megatree_runner.run()
    """

    def __init__(self, config: Config):
        """
        Initialize the MegatreeLoaderRunner with a configuration object.

        :param config: Configuration object containing tool and logger settings.
        :type config: Config
        """
        super().__init__(config)
        self.tool_name = 'megatree-loader'

    def set_infile(self, infile: str) -> None:
        """
        Set the input tree file in Newick format.

        :param infile: Path to the input file.
        :type infile: str
        """
        self.set_parameter('i', infile)

    def set_dbfile(self, dbfile: str) -> None:
        """
        Set the output database file.

        :param dbfile: Path to the output database file.
        :type dbfile: str
        """
        self.set_parameter('d', dbfile)

    def build_command(self) -> List[str]:
        """
        Build the megatree-loader command with all set parameters.

        :return: The complete megatree-loader command as a list of strings.
        :rtype: List[str]
        """
        command = super().build_command()

        # Handle verbosity based on logger level
        logger_level = self.logger.getEffectiveLevel()
        if logger_level <= logging.DEBUG:
            command.extend(['-v', '-v'])
        elif logger_level <= logging.INFO:
            command.append('-v')

        # Handle boolean flags (parameters without values)
        bool_params = ['h', 'm']
        for param in bool_params:
            if self.get_parameter(param) is not None:
                command.append(f'-{param}')

        return command

    def validate_parameters(self) -> None:
        """
        Validate that required parameters are set before running the command.

        :raises ValueError: If required parameters are missing.
        """
        if self.get_parameter('i') is None:
            raise ValueError("Input file (-i) is required for megatree-loader.")

    def run(self) -> int:
        """
        Run the megatree-loader command after validating parameters.

        :return: The return code of the megatree-loader command.
        :rtype: int
        :raises ValueError: If required parameters are missing.
        """
        self.validate_parameters()
        return super().run()
