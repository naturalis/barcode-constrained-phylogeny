from typing import List, Union
from bactria.tool_runner import ToolRunner
from bactria.config import Config


class BLASTNRunner(ToolRunner):
    """
    A subclass of ToolRunner specifically for running BLASTN (Nucleotide-Nucleotide BLAST).
    """

    def __init__(self, config: Config):
        """
        Initialize the BLASTNRunner with a configuration object.

        :param config: Configuration object containing tool and logger settings.
        :type config: Config
        """
        super().__init__(config)
        self.tool_name = 'blastn'
        self.set_task('megablast')  # Default task
        self.set_evalue(10)  # Default evalue
        self.set_outfmt(0)  # Default output format
        self.set_num_threads(1)  # Default number of threads
        self.set_max_target_seqs(500)  # Default max target sequences
        self.set_word_size(28)  # Default word size for megablast
        self.set_dust('20 64 1')  # Default DUST settings
        self.set_soft_masking(True)  # Default soft masking

    def set_query(self, query: str) -> None:
        """
        Set the input file name.

        :param query: Path to the input file.
        :type query: str
        """
        self.set_parameter('query', query)

    def set_db(self, db: str) -> None:
        """
        Set the BLAST database name.

        :param db: Name of the BLAST database.
        :type db: str
        """
        self.set_parameter('db', db)

    def set_out(self, out: str) -> None:
        """
        Set the output file name.

        :param out: Path to the output file.
        :type out: str
        """
        self.set_parameter('out', out)

    def set_evalue(self, evalue: float) -> None:
        """
        Set the expectation value (E) threshold for saving hits.

        :param evalue: E-value threshold.
        :type evalue: float
        """
        self.set_parameter('evalue', str(evalue))

    def set_outfmt(self, outfmt: Union[int, str]) -> None:
        """
        Set the alignment view options.

        :param outfmt: Output format (0-18 or custom specification).
        :type outfmt: Union[int, str]
        """
        self.set_parameter('outfmt', str(outfmt))

    def set_num_threads(self, num_threads: int) -> None:
        """
        Set the number of threads (CPUs) to use in the BLAST search.

        :param num_threads: Number of threads to use.
        :type num_threads: int
        """
        self.set_parameter('num_threads', str(num_threads))

    def set_max_target_seqs(self, max_target_seqs: int) -> None:
        """
        Set the maximum number of aligned sequences to keep.

        :param max_target_seqs: Maximum number of aligned sequences.
        :type max_target_seqs: int
        """
        self.set_parameter('max_target_seqs', str(max_target_seqs))

    def set_task(self, task: str) -> None:
        """
        Set the task to execute.

        :param task: Task name (e.g., 'blastn', 'megablast', 'dc-megablast').
        :type task: str
        """
        self.set_parameter('task', task)

    def set_word_size(self, word_size: int) -> None:
        """
        Set the word size for wordfinder algorithm.

        :param word_size: Word size (>=4).
        :type word_size: int
        """
        self.set_parameter('word_size', str(word_size))

    def set_dust(self, dust: str) -> None:
        """
        Set the DUST filtering options.

        :param dust: DUST options ('yes', 'level window linker', or 'no').
        :type dust: str
        """
        self.set_parameter('dust', dust)

    def set_soft_masking(self, soft_masking: bool) -> None:
        """
        Set whether to apply filtering locations as soft masks.

        :param soft_masking: Whether to use soft masking.
        :type soft_masking: bool
        """
        self.set_parameter('soft_masking', str(soft_masking).lower())

    def set_perc_identity(self, perc_identity: float) -> None:
        """
        Set the percent identity threshold.

        :param perc_identity: Percent identity (0-100).
        :type perc_identity: float
        """
        self.set_parameter('perc_identity', str(perc_identity))

    def set_culling_limit(self, culling_limit: int) -> None:
        """
        Set the culling limit.

        :param culling_limit: Culling limit (>=0).
        :type culling_limit: int
        """
        self.set_parameter('culling_limit', str(culling_limit))

    def set_best_hit_overhang(self, best_hit_overhang: float) -> None:
        """
        Set the Best Hit algorithm overhang value.

        :param best_hit_overhang: Best hit overhang value (>0 and <0.5).
        :type best_hit_overhang: float
        """
        self.set_parameter('best_hit_overhang', str(best_hit_overhang))

    def set_best_hit_score_edge(self, best_hit_score_edge: float) -> None:
        """
        Set the Best Hit algorithm score edge value.

        :param best_hit_score_edge: Best hit score edge value (>0 and <0.5).
        :type best_hit_score_edge: float
        """
        self.set_parameter('best_hit_score_edge', str(best_hit_score_edge))

    def set_taxids(self, taxids: List[str]) -> None:
        """
        Set the taxonomy IDs to restrict the search.

        :param taxids: List of taxonomy IDs.
        :type taxids: List[str]
        """
        self.set_parameter('taxids', ','.join(taxids))

    def set_negative_taxids(self, negative_taxids: List[str]) -> None:
        """
        Set the taxonomy IDs to exclude from the search.

        :param negative_taxids: List of taxonomy IDs to exclude.
        :type negative_taxids: List[str]
        """
        self.set_parameter('negative_taxids', ','.join(negative_taxids))

    def set_entrez_query(self, entrez_query: str) -> None:
        """
        Set the Entrez query to restrict the search.

        :param entrez_query: Entrez query string.
        :type entrez_query: str
        """
        self.set_parameter('entrez_query', entrez_query)

    def set_remote(self, remote: bool = True) -> None:
        """
        Set whether to execute the search remotely.

        :param remote: Whether to run the search remotely.
        :type remote: bool
        """
        if remote:
            self.set_parameter('remote', '')
        else:
            self.parameters.pop('remote', None)

    def build_command(self) -> List[str]:
        """
        Build the BLASTN command with all set parameters.

        :return: The complete BLASTN command as a list of strings.
        :rtype: List[str]
        """
        command = super().build_command()

        # Handle boolean flags (parameters without values)
        bool_params = ['remote', 'parse_deflines', 'show_gis', 'ungapped', 'use_index', 'lcase_masking']
        for param in bool_params:
            if self.get_parameter(param) is not None:
                command.append(f'-{param}')

        return command

# Example usage:
# config = Config()
# config.load_config('path/to/config.yaml')
# blastn_runner = BLASTNRunner(config)
# blastn_runner.set_query('input.fasta')
# blastn_runner.set_db('nr')
# blastn_runner.set_out('output.txt')
# blastn_runner.set_evalue(1e-5)
# blastn_runner.set_outfmt(6)
# blastn_runner.set_num_threads(4)
# blastn_runner.set_remote(True)
# return_code = blastn_runner.run()