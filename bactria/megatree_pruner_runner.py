from typing import List
from bactria.tool_runner import ToolRunner
from bactria.config import Config
import logging


class MegatreePrunerRunner(ToolRunner):
    """
    A subclass of ToolRunner specifically for running megatree-pruner.
    """

    def __init__(self, config: Config):
        """
        Initialize the MegatreePrunerRunner with a configuration object.

        :param config: Configuration object containing tool and logger settings.
        :type config: Config
        """
        super().__init__(config)
        self.tool_name = 'megatree-pruner'

    def set_dbfile(self, dbfile: str) -> None:
        """
        Set the database file produced by megatree-*-loader scripts.

        :param dbfile: Path to the database file.
        :type dbfile: str
        """
        self.set_parameter('d', dbfile)

    def set_infile(self, infile: str) -> None:
        """
        Set the input file containing a list of taxon names.

        :param infile: Path to the input file.
        :type infile: str
        """
        self.set_parameter('i', infile)

    def set_list(self, taxa_list: List[str]) -> None:
        """
        Set the list of taxon names to be retained.

        :param taxa_list: List of taxon names.
        :type taxa_list: List[str]
        """
        self.set_parameter('l', ','.join(taxa_list))

    def set_tabular(self, tabular: bool = True) -> None:
        """
        Set the option to produce a tab-separated table instead of a Newick-formatted tree.

        :param tabular: Whether to produce tabular output.
        :type tabular: bool
        """
        if tabular:
            self.set_parameter('t', '')

    def set_relabel(self, relabel: bool = True) -> None:
        """
        Set the option to relabel internal nodes in the output.

        :param relabel: Whether to relabel internal nodes.
        :type relabel: bool
        """
        if relabel:
            self.set_parameter('r', '')

    def build_command(self) -> List[str]:
        """
        Build the megatree-pruner command with all set parameters.

        :return: The complete megatree-pruner command as a list of strings.
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
        bool_params = ['t', 'r', 'h', 'm']
        for param in bool_params:
            if self.get_parameter(param) is not None:
                command.append(f'-{param}')

        return command

    def validate_parameters(self) -> None:
        """
        Validate that required parameters are set before running the command.

        :raises ValueError: If required parameters are missing or incompatible options are set.
        """
        if self.get_parameter('d') is None:
            raise ValueError("Database file (-d) is required for megatree-pruner.")

        if self.get_parameter('i') is None and self.get_parameter('l') is None:
            raise ValueError("Either input file (-i) or taxon list (-l) must be provided for megatree-pruner.")

        if self.get_parameter('i') is not None and self.get_parameter('l') is not None:
            raise ValueError("Only one of input file (-i) or taxon list (-l) can be provided for megatree-pruner.")

    def run(self) -> int:
        """
        Run the megatree-pruner command after validating parameters.

        :return: The return code of the megatree-pruner command.
        :rtype: int
        :raises ValueError: If required parameters are missing or incompatible options are set.
        """
        self.validate_parameters()
        return super().run()

# Example usage:
# config = Config()
# config.load_config('path/to/config.yaml')
# megatree_pruner = MegatreePrunerRunner(config)
# megatree_pruner.set_dbfile('tree_database.sqlite')
# megatree_pruner.set_infile('taxa_list.txt')
# megatree_pruner.set_tabular(True)
# megatree_pruner.set_relabel(True)
# return_code = megatree_pruner.run()