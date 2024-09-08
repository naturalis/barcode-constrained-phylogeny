from typing import Optional, List
from bactria.tool_runner import ToolRunner
from bactria.config import Config


class Blastdbcommand(ToolRunner):
    """
    A subclass of ToolRunner specifically for running blastdbcmd.

    Examples:
        >>> config = Config()
        >>> config.load_config('path/to/config.yaml')
        >>> blastdbcmd_runner = Blastdbcommand(config)
        >>> blastdbcmd_runner.set_db('nr')
        >>> blastdbcmd_runner.set_entry(['NP_000508.1', 'NP_001018081.2'])
        >>> blastdbcmd_runner.set_outfmt('%f')
        >>> return_code = blastdbcmd_runner.run()
    """

    def __init__(self, config: Config):
        """
        Initialize the BlastdbcommandRunner with a configuration object.

        :param config: Configuration object containing tool and logger settings.
        :type config: Config
        """
        super().__init__(config)
        self.tool_name = 'blastdbcmd'
        self.set_db('nr')  # Default database
        self.set_dbtype('guess')  # Default dbtype

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

        :param molecule_type: Molecule type ('guess', 'nucl', or 'prot').
        :type molecule_type: str
        """
        self.set_parameter('dbtype', molecule_type)

    def set_entry(self, identifiers: List[str]) -> None:
        """
        Set the comma-delimited search string(s) of sequence identifiers.

        :param identifiers: List of sequence identifiers.
        :type identifiers: List[str]
        """
        self.set_parameter('entry', ','.join(identifiers))

    def set_entry_batch(self, file_path: str) -> None:
        """
        Set the input file for batch processing.

        :param file_path: Path to the input file.
        :type file_path: str
        """
        self.set_parameter('entry_batch', file_path)

    def set_ipg(self, ipg: int) -> None:
        """
        Set the IPG to retrieve.

        :param ipg: IPG number.
        :type ipg: int
        """
        self.set_parameter('ipg', str(ipg))

    def set_ipg_batch(self, file_path: str) -> None:
        """
        Set the input file for batch processing IPGs.

        :param file_path: Path to the input file.
        :type file_path: str
        """
        self.set_parameter('ipg_batch', file_path)

    def set_info(self) -> None:
        """
        Set the flag to print BLAST database information.
        """
        self.set_parameter('info', '')

    def set_metadata(self) -> None:
        """
        Set the flag to generate BLAST database metadata.
        """
        self.set_parameter('metadata', '')

    def set_metadata_output_prefix(self, prefix: str) -> None:
        """
        Set the path prefix for location of database files in metadata.

        :param prefix: Path prefix.
        :type prefix: str
        """
        self.set_parameter('metadata_output_prefix', prefix)

    def set_tax_info(self) -> None:
        """
        Set the flag to print taxonomic information.
        """
        self.set_parameter('tax_info', '')

    def set_range(self, start: int, stop: Optional[int] = None) -> None:
        """
        Set the range of sequence to extract in 1-based offsets.

        :param start: Start position.
        :type start: int
        :param stop: Stop position (optional).
        :type stop: Optional[int]
        """
        range_str = f"{start}-{stop if stop is not None else ''}"
        self.set_parameter('range', range_str)

    def set_strand(self, strand: str) -> None:
        """
        Set the strand of nucleotide sequence to extract.

        :param strand: Strand ('minus' or 'plus').
        :type strand: str
        """
        self.set_parameter('strand', strand)

    def set_mask_sequence_with(self, algo_id: str) -> None:
        """
        Set the algorithm ID for producing lower-case masked FASTA.

        :param algo_id: Algorithm ID.
        :type algo_id: str
        """
        self.set_parameter('mask_sequence_with', algo_id)

    def set_taxids(self, tax_ids: List[str]) -> None:
        """
        Set the comma-delimited taxonomy identifiers.

        :param tax_ids: List of taxonomy identifiers.
        :type tax_ids: List[str]
        """
        self.set_parameter('taxids', ','.join(tax_ids))

    def set_taxidlist(self, file_path: str) -> None:
        """
        Set the input file for taxonomy identifiers.

        :param file_path: Path to the input file.
        :type file_path: str
        """
        self.set_parameter('taxidlist', file_path)

    def set_no_taxid_expansion(self) -> None:
        """
        Set the flag to not expand the taxonomy IDs to their descendant taxonomy IDs.
        """
        self.set_parameter('no_taxid_expansion', '')

    def set_out(self, file_path: str) -> None:
        """
        Set the output file name.

        :param file_path: Path to the output file.
        :type file_path: str
        """
        self.set_parameter('out', file_path)

    def set_outfmt(self, format_string: str) -> None:
        """
        Set the output format.

        :param format_string: Output format string.
        :type format_string: str
        """
        self.set_parameter('outfmt', format_string)

    def set_target_only(self) -> None:
        """
        Set the flag for definition line to contain target entry only.
        """
        self.set_parameter('target_only', '')

    def set_get_dups(self) -> None:
        """
        Set the flag to retrieve duplicate accessions.
        """
        self.set_parameter('get_dups', '')

    def set_line_length(self, length: int) -> None:
        """
        Set the line length for output.

        :param length: Line length.
        :type length: int
        """
        self.set_parameter('line_length', str(length))

    def set_ctrl_a(self) -> None:
        """
        Set the flag to use Ctrl-A as the non-redundant defline separator.
        """
        self.set_parameter('ctrl_a', '')

    def set_show_blastdb_search_path(self) -> None:
        """
        Set the flag to display the default BLAST database search paths.
        """
        self.set_parameter('show_blastdb_search_path', '')

    def set_list(self, directory: str) -> None:
        """
        Set the directory to list BLAST databases.

        :param directory: Directory path.
        :type directory: str
        """
        self.set_parameter('list', directory)

    def set_remove_redundant_dbs(self) -> None:
        """
        Set the flag to remove redundant databases.
        """
        self.set_parameter('remove_redundant_dbs', '')

    def set_recursive(self) -> None:
        """
        Set the flag to recursively traverse the directory structure.
        """
        self.set_parameter('recursive', '')

    def set_list_outfmt(self, format_string: str) -> None:
        """
        Set the output format for the list option.

        :param format_string: Output format string.
        :type format_string: str
        """
        self.set_parameter('list_outfmt', format_string)

    def set_exact_length(self) -> None:
        """
        Set the flag to get exact length for db info.
        """
        self.set_parameter('exact_length', '')

    def set_long_seqids(self) -> None:
        """
        Set the flag to use long seq id for fasta deflines.
        """
        self.set_parameter('long_seqids', '')

    def set_logfile(self, file_path: str) -> None:
        """
        Set the file to which the program log should be redirected.

        :param file_path: Path to the log file.
        :type file_path: str
        """
        self.set_parameter('logfile', file_path)

    def validate_parameters(self) -> None:
        """
        Validate that required parameters are set and incompatible options are not used together.

        :raises ValueError: If required parameters are missing or incompatible options are set.
        """
        # This method should implement validation logic based on the manual page.
        # For brevity, only a few checks are implemented here.
        if self.get_parameter('entry') and self.get_parameter('entry_batch'):
            raise ValueError("'entry' and 'entry_batch' options are incompatible.")

        if self.get_parameter('info') and self.get_parameter('outfmt'):
            raise ValueError("'info' and 'outfmt' options are incompatible.")

    def run(self) -> int:
        """
        Run the blastdbcmd command after validating parameters.

        :return: The return code of the blastdbcmd command.
        :rtype: int
        :raises ValueError: If required parameters are missing or incompatible options are set.
        """
        self.validate_parameters()
        return super().run()
