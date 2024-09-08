from typing import List
from .tool_runner import ToolRunner
from bactria.config import Config


class Hmmalign(ToolRunner):
    """
    A subclass of ToolRunner specifically for running hmmalign.

    Examples:
        >>> config = Config()
        >>> config.load_config('path/to/config.yaml')
        >>> hmmalign_runner = Hmmalign(config)
        >>> hmmalign_runner.set_hmmfile('model.hmm')
        >>> hmmalign_runner.set_seqfile('sequences.fasta')
        >>> hmmalign_runner.set_output('alignment.sto')
        >>> hmmalign_runner.set_outformat('Stockholm')
        >>> return_code = hmmalign_runner.run()
    """

    def __init__(self, config: Config):
        """
        Initialize the HmmalignRunner with a configuration object.

        :param config: Configuration object containing tool and logger settings.
        :type config: Config
        """
        super().__init__(config)
        self.tool_name = 'hmmalign'

    def set_hmmfile(self, hmmfile: str) -> None:
        """
        Set the HMM file.

        :param hmmfile: Path to the HMM file.
        :type hmmfile: str
        """
        self.set_parameter('hmmfile', hmmfile)

    def set_seqfile(self, seqfile: str) -> None:
        """
        Set the sequence file.

        :param seqfile: Path to the sequence file.
        :type seqfile: str
        """
        self.set_parameter('seqfile', seqfile)

    def set_output(self, output_file: str) -> None:
        """
        Set the output file for the alignment.

        :param output_file: Path to the output file.
        :type output_file: str
        """
        self.set_parameter('o', output_file)

    def set_mapali(self, mapali_file: str) -> None:
        """
        Set the file containing the alignment that the HMM came from.

        :param mapali_file: Path to the mapali file.
        :type mapali_file: str
        """
        self.set_parameter('mapali', mapali_file)

    def set_trim(self) -> None:
        """
        Set the flag to trim terminal tails of nonaligned residues from alignment.
        """
        self.set_parameter('trim', '')

    def set_amino(self) -> None:
        """
        Set the flag to assert that seqfile and hmmfile are both protein.
        """
        self.set_parameter('amino', '')

    def set_dna(self) -> None:
        """
        Set the flag to assert that seqfile and hmmfile are both DNA.
        """
        self.set_parameter('dna', '')

    def set_rna(self) -> None:
        """
        Set the flag to assert that seqfile and hmmfile are both RNA.
        """
        self.set_parameter('rna', '')

    def set_informat(self, format: str) -> None:
        """
        Set the input format for the sequence file.

        :param format: Input format (e.g., 'FASTA', 'EMBL', 'GenBank', 'UniProt').
        :type format: str
        """
        self.set_parameter('informat', format)

    def set_outformat(self, format: str) -> None:
        """
        Set the output format for the alignment.

        :param format: Output format (e.g., 'Stockholm', 'Pfam', 'A2M', 'PSIBLAST').
        :type format: str
        """
        self.set_parameter('outformat', format)

    def validate_parameters(self) -> None:
        """
        Validate that required parameters are set and incompatible options are not used together.

        :raises ValueError: If required parameters are missing or incompatible options are set.
        """
        if not self.get_parameter('hmmfile'):
            raise ValueError("HMM file is required for hmmalign.")
        if not self.get_parameter('seqfile'):
            raise ValueError("Sequence file is required for hmmalign.")

        # Check for mutually exclusive options
        molecule_types = ['amino', 'dna', 'rna']
        set_types = sum(1 for t in molecule_types if self.get_parameter(t) is not None)
        if set_types > 1:
            raise ValueError("Only one of --amino, --dna, or --rna can be set.")

    def build_command(self) -> List[str]:
        """
        Build the hmmalign command with all set parameters.

        :return: The complete hmmalign command as a list of strings.
        :rtype: List[str]
        """
        command = super().build_command()

        # Add positional arguments at the end
        if self.get_parameter('hmmfile'):
            command.append(self.get_parameter('hmmfile'))
        if self.get_parameter('seqfile'):
            command.append(self.get_parameter('seqfile'))

        return command

    def run(self) -> int:
        """
        Run the hmmalign command after validating parameters.

        :return: The return code of the hmmalign command.
        :rtype: int
        :raises ValueError: If required parameters are missing or incompatible options are set.
        """
        self.validate_parameters()
        return super().run()
