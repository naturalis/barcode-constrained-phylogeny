from typing import Optional, List, Union
from .tool_runner import ToolRunner
from bactria.config import Config


class RaxmlNg(ToolRunner):
    """
    A subclass of ToolRunner specifically for running RAxML-NG (Randomized Axelerated Maximum Likelihood - Next Generation).

    Examples:
        >>> config = Config()
        >>> config.load_config('path/to/config.yaml')
        >>> raxml_runner = RaxmlNg(config)
        >>> raxml_runner.set_msa('alignment.fa')
        >>> raxml_runner.set_model('GTR+G')
        >>> raxml_runner.set_prefix('test_run')
        >>> raxml_runner.set_threads(4)
        >>> raxml_runner.set_all()
        >>> return_code = raxml_runner.run()
    """

    def __init__(self, config: Config):
        """
        Initialize the RAXMLNGRunner with a configuration object.

        :param config: Configuration object containing tool and logger settings.
        :type config: Config
        """
        super().__init__(config)
        self.tool_name = 'raxml-ng'
        self.set_threads(6)  # Default number of threads
        self.set_workers(1)  # Default number of workers
        self.set_seed(None)  # Default seed (current time)
        self.set_blmin(1e-6)  # Default minimum branch length
        self.set_blmax(100)  # Default maximum branch length
        self.set_lh_epsilon(0.1)  # Default log-likelihood epsilon

    def set_msa(self, msa: str) -> None:
        """
        Set the alignment file.

        :param msa: Path to the alignment file.
        :type msa: str
        """
        self.set_parameter('msa', msa)

    def set_model(self, model: str) -> None:
        """
        Set the model specification or partition file.

        :param model: Model specification or path to partition file.
        :type model: str
        """
        self.set_parameter('model', model)

    def set_prefix(self, prefix: str) -> None:
        """
        Set the prefix for output files.

        :param prefix: Prefix for output files.
        :type prefix: str
        """
        self.set_parameter('prefix', prefix)

    def set_threads(self, threads: int) -> None:
        """
        Set the number of parallel threads to use.

        :param threads: Number of threads.
        :type threads: int
        """
        self.set_parameter('threads', str(threads))

    def set_workers(self, workers: int) -> None:
        """
        Set the number of tree searches to run in parallel.

        :param workers: Number of workers.
        :type workers: int
        """
        self.set_parameter('workers', str(workers))

    def set_seed(self, seed: Optional[int]) -> None:
        """
        Set the seed for pseudo-random number generator.

        :param seed: Seed value (None for current time).
        :type seed: Optional[int]
        """
        if seed is not None:
            self.set_parameter('seed', str(seed))

    def set_tree(self, tree: str) -> None:
        """
        Set the starting tree option.

        :param tree: Starting tree option (e.g., 'rand{10}', 'pars{10}', or file path).
        :type tree: str
        """
        self.set_parameter('tree', tree)

    def set_blmin(self, blmin: float) -> None:
        """
        Set the minimum branch length.

        :param blmin: Minimum branch length.
        :type blmin: float
        """
        self.set_parameter('blmin', str(blmin))

    def set_blmax(self, blmax: float) -> None:
        """
        Set the maximum branch length.

        :param blmax: Maximum branch length.
        :type blmax: float
        """
        self.set_parameter('blmax', str(blmax))

    def set_lh_epsilon(self, lh_epsilon: float) -> None:
        """
        Set the log-likelihood epsilon for optimization/tree search.

        :param lh_epsilon: Log-likelihood epsilon.
        :type lh_epsilon: float
        """
        self.set_parameter('lh-epsilon', str(lh_epsilon))

    def set_brlen(self, brlen: str) -> None:
        """
        Set the branch length linkage between partitions.

        :param brlen: Branch length linkage ('linked', 'scaled', or 'unlinked').
        :type brlen: str
        """
        self.set_parameter('brlen', brlen)

    def set_spr_radius(self, spr_radius: Union[int, str]) -> None:
        """
        Set the SPR re-insertion radius for fast iterations.

        :param spr_radius: SPR radius (int or 'AUTO').
        :type spr_radius: Union[int, str]
        """
        self.set_parameter('spr-radius', str(spr_radius))

    def set_bs_trees(self, bs_trees: Union[int, str]) -> None:
        """
        Set the number of bootstrap replicates or autoMRE option.

        :param bs_trees: Number of bootstrap replicates or 'autoMRE{N}'.
        :type bs_trees: Union[int, str]
        """
        self.set_parameter('bs-trees', str(bs_trees))

    def set_all(self) -> None:
        """
        Set the all-in-one analysis option (ML search + bootstrapping).
        """
        self.set_parameter('all', '')

    def set_search(self) -> None:
        """
        Set the ML tree search option.
        """
        self.set_parameter('search', '')

    def set_evaluate(self) -> None:
        """
        Set the evaluate option (likelihood of a tree with model+brlen optimization).
        """
        self.set_parameter('evaluate', '')

    def build_command(self) -> List[str]:
        """
        Build the RAxML-NG command with all set parameters.

        :return: The complete RAxML-NG command as a list of strings.
        :rtype: List[str]
        """
        command = super().build_command()

        # Handle boolean flags (parameters without values)
        bool_params = ['all', 'search', 'evaluate', 'bootstrap', 'support', 'redo', 'nofiles']
        for param in bool_params:
            if self.get_parameter(param) is not None:
                command.append(f'--{param}')

        return command
