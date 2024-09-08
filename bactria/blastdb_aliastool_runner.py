from typing import Optional, List
from bactria.tool_runner import ToolRunner
from bactria.config import Config


class BlastdbAliastool(ToolRunner):
    """
    A subclass of ToolRunner specifically for running blastdb_aliastool.

    Examples:
        >>> config = Config()
        >>> config.load_config('path/to/config.yaml')
        >>> aliastool_runner = BlastdbAliastoolRunner(config)
        >>> aliastool_runner.set_db('my_blast_db')
        >>> aliastool_runner.set_out('my_alias_db')
        >>> aliastool_runner.set_gilist('my_gi_list.txt')
        >>> return_code = aliastool_runner.run()
    """

    def __init__(self, config: Config):
        """
        Initialize the BlastdbAliastoolRunner with a configuration object.

        :param config: Configuration object containing tool and logger settings.
        :type config: Config
        """
        super().__init__(config)
        self.tool_name = 'blastdb_aliastool'

    def set_gi_file_in(self, file_path: str) -> None:
        """
        Set the input text file containing GIs to convert.

        :param file_path: Path to the input file.
        :type file_path: str
        """
        self.set_parameter('gi_file_in', file_path)

    def set_gi_file_out(self, file_path: str) -> None:
        """
        Set the output file name for the converted GI file.

        :param file_path: Path to the output file.
        :type file_path: str
        """
        self.set_parameter('gi_file_out', file_path)

    def set_db(self, dbname: str) -> None:
        """
        Set the BLAST database name.

        :param dbname: Name of the BLAST database.
        :type dbname: str
        """
        self.set_parameter('db', dbname)

    def set_dbtype(self, molecule_type: str) -> None:
        """
        Set the molecule type stored in the BLAST database.

        :param molecule_type: Molecule type ('nucl' or 'prot').
        :type molecule_type: str
        """
        self.set_parameter('dbtype', molecule_type)

    def set_title(self, title: str) -> None:
        """
        Set the title for the BLAST database.

        :param title: Title for the database.
        :type title: str
        """
        self.set_parameter('title', title)

    def set_gilist(self, file_path: str) -> None:
        """
        Set the text or binary gi file to restrict the BLAST database.

        :param file_path: Path to the gi file.
        :type file_path: str
        """
        self.set_parameter('gilist', file_path)

    def set_seqidlist(self, file_path: str) -> None:
        """
        Set the text sequence id or accession file to restrict the BLAST database.

        :param file_path: Path to the sequence id file.
        :type file_path: str
        """
        self.set_parameter('seqidlist', file_path)

    def set_taxidlist(self, file_path: str) -> None:
        """
        Set the text taxonomy id file to restrict the BLAST database.

        :param file_path: Path to the taxonomy id file.
        :type file_path: str
        """
        self.set_parameter('taxidlist', file_path)

    def set_oid_masks(self, mask_value: int) -> None:
        """
        Set the oid masks for creating alias db.

        :param mask_value: Integer value for oid masks.
        :type mask_value: int
        """
        self.set_parameter('oid_masks', str(mask_value))

    def set_out(self, database_name: str) -> None:
        """
        Set the name of the BLAST database alias to be created.

        :param database_name: Name of the database alias.
        :type database_name: str
        """
        self.set_parameter('out', database_name)

    def set_dblist(self, database_names: List[str]) -> None:
        """
        Set a list of BLAST database names to aggregate.

        :param database_names: List of database names.
        :type database_names: List[str]
        """
        self.set_parameter('dblist', ' '.join(database_names))

    def set_dblist_file(self, file_path: str) -> None:
        """
        Set a file containing a list of BLAST database names to aggregate.

        :param file_path: Path to the file containing database names.
        :type file_path: str
        """
        self.set_parameter('dblist_file', file_path)

    def set_vdblist(self, vdb_names: List[str]) -> None:
        """
        Set a list of VDB names to aggregate.

        :param vdb_names: List of VDB names.
        :type vdb_names: List[str]
        """
        self.set_parameter('vdblist', ' '.join(vdb_names))

    def set_vdblist_file(self, file_path: str) -> None:
        """
        Set a file containing a list of VDB names to aggregate.

        :param file_path: Path to the file containing VDB names.
        :type file_path: str
        """
        self.set_parameter('vdblist_file', file_path)

    def set_num_volumes(self, num: int) -> None:
        """
        Set the number of volumes to aggregate.

        :param num: Number of volumes.
        :type num: int
        """
        self.set_parameter('num_volumes', str(num))

    def set_logfile(self, file_path: str) -> None:
        """
        Set the file to which the program log should be redirected.

        :param file_path: Path to the log file.
        :type file_path: str
        """
        self.set_parameter('logfile', file_path)

    def set_version(self) -> None:
        """
        Set the version flag to print the version number.
        """
        self.set_parameter('version', '')

    def validate_parameters(self) -> None:
        """
        Validate that required parameters are set and incompatible options are not used together.

        :raises ValueError: If required parameters are missing or incompatible options are set.
        """
        # This method should implement validation logic based on the manual page.
        # For brevity, only a few checks are implemented here.
        if self.get_parameter('out') and self.get_parameter('gi_file_in'):
            raise ValueError("'out' and 'gi_file_in' options are incompatible.")

        if self.get_parameter('dblist') and not (
                self.get_parameter('out') and self.get_parameter('dbtype') and self.get_parameter('title')):
            raise ValueError("'dblist' requires 'out', 'dbtype', and 'title' to be set.")

    def run(self) -> int:
        """
        Run the blastdb_aliastool command after validating parameters.

        :return: The return code of the blastdb_aliastool command.
        :rtype: int
        :raises ValueError: If required parameters are missing or incompatible options are set.
        """
        self.validate_parameters()
        return super().run()
