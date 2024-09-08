from typing import List
from bactria.tool_runner import ToolRunner
from bactria.config import Config

class MakeblastdbRunner(ToolRunner):
    """
    A subclass of ToolRunner specifically for running makeblastdb.
    """

    def __init__(self, config: Config):
        """
        Initialize the MakeblastdbRunner with a configuration object.

        :param config: Configuration object containing tool and logger settings.
        :type config: Config
        """
        super().__init__(config)
        self.tool_name = 'makeblastdb'
        self.set_input_type('fasta')  # Default input type
        self.set_blastdb_version(5)  # Default BLAST database version
        self.set_max_file_sz('3GB')  # Default maximum file size

    def set_in(self, input_file: str) -> None:
        """
        Set the input file/database name.

        :param input_file: Path to the input file.
        :type input_file: str
        """
        self.set_parameter('in', input_file)

    def set_input_type(self, input_type: str) -> None:
        """
        Set the type of the data specified in input_file.

        :param input_type: Type of input data ('asn1_bin', 'asn1_txt', 'blastdb', 'fasta').
        :type input_type: str
        """
        self.set_parameter('input_type', input_type)

    def set_dbtype(self, molecule_type: str) -> None:
        """
        Set the molecule type of target db.

        :param molecule_type: Molecule type ('nucl' or 'prot').
        :type molecule_type: str
        """
        self.set_parameter('dbtype', molecule_type)

    def set_title(self, title: str) -> None:
        """
        Set the title for BLAST database.

        :param title: Title for the database.
        :type title: str
        """
        self.set_parameter('title', title)

    def set_parse_seqids(self) -> None:
        """
        Set the flag to parse seqid for FASTA input.
        """
        self.set_parameter('parse_seqids', '')

    def set_hash_index(self) -> None:
        """
        Set the flag to create index of sequence hash values.
        """
        self.set_parameter('hash_index', '')

    def set_mask_data(self, mask_data_files: List[str]) -> None:
        """
        Set the list of input files containing masking data.

        :param mask_data_files: List of paths to masking data files.
        :type mask_data_files: List[str]
        """
        self.set_parameter('mask_data', ','.join(mask_data_files))

    def set_mask_id(self, mask_algo_ids: List[str]) -> None:
        """
        Set the list of strings to uniquely identify the masking algorithm.

        :param mask_algo_ids: List of masking algorithm identifiers.
        :type mask_algo_ids: List[str]
        """
        self.set_parameter('mask_id', ','.join(mask_algo_ids))

    def set_mask_desc(self, mask_algo_descriptions: List[str]) -> None:
        """
        Set the list of strings to describe the masking algorithm details.

        :param mask_algo_descriptions: List of masking algorithm descriptions.
        :type mask_algo_descriptions: List[str]
        """
        self.set_parameter('mask_desc', ','.join(mask_algo_descriptions))

    def set_gi_mask(self) -> None:
        """
        Set the flag to create GI indexed masking data.
        """
        self.set_parameter('gi_mask', '')

    def set_gi_mask_name(self, gi_based_mask_names: List[str]) -> None:
        """
        Set the list of masking data output files.

        :param gi_based_mask_names: List of masking data output file names.
        :type gi_based_mask_names: List[str]
        """
        self.set_parameter('gi_mask_name', ','.join(gi_based_mask_names))

    def set_out(self, database_name: str) -> None:
        """
        Set the name of BLAST database to be created.

        :param database_name: Name of the BLAST database.
        :type database_name: str
        """
        self.set_parameter('out', database_name)

    def set_blastdb_version(self, version: int) -> None:
        """
        Set the version of BLAST database to be created.

        :param version: BLAST database version (4 or 5).
        :type version: int
        """
        self.set_parameter('blastdb_version', str(version))

    def set_max_file_sz(self, number_of_bytes: str) -> None:
        """
        Set the maximum file size for BLAST database files.

        :param number_of_bytes: Maximum file size (e.g., '3GB').
        :type number_of_bytes: str
        """
        self.set_parameter('max_file_sz', number_of_bytes)

    def set_metadata_output_prefix(self, prefix: str) -> None:
        """
        Set the path prefix for location of database files in metadata.

        :param prefix: Path prefix for metadata.
        :type prefix: str
        """
        self.set_parameter('metadata_output_prefix', prefix)

    def set_logfile(self, file_name: str) -> None:
        """
        Set the file to which the program log should be redirected.

        :param file_name: Path to the log file.
        :type file_name: str
        """
        self.set_parameter('logfile', file_name)

    def set_taxid(self, taxid: int) -> None:
        """
        Set the taxonomy ID to assign to all sequences.

        :param taxid: Taxonomy ID.
        :type taxid: int
        """
        self.set_parameter('taxid', str(taxid))

    def set_taxid_map(self, taxid_map_file: str) -> None:
        """
        Set the text file mapping sequence IDs to taxonomy IDs.

        :param taxid_map_file: Path to the taxid map file.
        :type taxid_map_file: str
        """
        self.set_parameter('taxid_map', taxid_map_file)

    def set_oid_masks(self, oid_masks: int) -> None:
        """
        Set the OID masks.

        :param oid_masks: OID masks value.
        :type oid_masks: int
        """
        self.set_parameter('oid_masks', str(oid_masks))

    def validate_parameters(self) -> None:
        """
        Validate that required parameters are set and incompatible options are not used together.

        :raises ValueError: If required parameters are missing or incompatible options are set.
        """
        if not self.get_parameter('dbtype'):
            raise ValueError("'dbtype' is a required parameter for makeblastdb.")

        if self.get_parameter('mask_id') and not self.get_parameter('mask_data'):
            raise ValueError("'mask_id' requires 'mask_data' to be set.")

        if self.get_parameter('mask_desc') and not self.get_parameter('mask_id'):
            raise ValueError("'mask_desc' requires 'mask_id' to be set.")

        if self.get_parameter('gi_mask') and self.get_parameter('mask_id'):
            raise ValueError("'gi_mask' and 'mask_id' are incompatible options.")

        if self.get_parameter('taxid') and self.get_parameter('taxid_map'):
            raise ValueError("'taxid' and 'taxid_map' are incompatible options.")

    def run(self) -> int:
        """
        Run the makeblastdb command after validating parameters.

        :return: The return code of the makeblastdb command.
        :rtype: int
        :raises ValueError: If required parameters are missing or incompatible options are set.
        """
        self.validate_parameters()
        return super().run()

# Example usage:
# config = Config()
# config.load_config('path/to/config.yaml')
# makeblastdb_runner = MakeblastdbRunner(config)
# makeblastdb_runner.set_in('input_sequences.fasta')
# makeblastdb_runner.set_dbtype('prot')
# makeblastdb_runner.set_out('my_blast_db')
# makeblastdb_runner.set_parse_seqids()
# return_code = makeblastdb_runner.run()