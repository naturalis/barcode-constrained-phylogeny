import os
import yaml
import copy


class Config:

    def __init__(self):
        self.config_data = None
        self.config_path = None
        self.initialized = False

    def load_config(self, config_path):
        """
        Load a configuration file.
        :param config_path:
        :return:
        """

        self.config_path = config_path
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Config file not found: {config_path}")

        with open(config_path, 'r') as config_file:
            try:
                self.config_data = yaml.safe_load(config_file)
            except yaml.YAMLError as e:
                raise yaml.YAMLError(f"Error parsing config file: {e}")

        self._process_relative_paths(os.path.dirname(os.path.abspath(config_path)))
        self.initialized = True

    def _process_relative_paths(self, config_dir):
        """
        Process relative paths in the configuration file.
        :param config_dir:
        :return:
        """
        for key, value in self.config_data.items():
            if isinstance(value, str) and not os.path.isabs(value) and key.endswith('_file'):
                self.config_data[key] = os.path.join(config_dir, value)

    def get(self, key, default=None):
        """
        Get a configuration value by key.
        :param key:
        :param default:
        :return:
        """
        if not self.initialized:
            raise RuntimeError("Configuration not loaded. Call load_config first.")
        return self.config_data.get(key, default)

    def set(self, key, value):
        """
        Set a configuration value by key.
        :param key:
        :param value:
        :return:
        """
        if not self.initialized:
            raise RuntimeError("Configuration not loaded. Call load_config first.")
        self.config_data[key] = value

    def detach(self):
        """
        Create a deep copy of the configuration object.
        :return:
        """
        if not self.initialized:
            raise RuntimeError("Configuration not loaded. Call load_config first.")
        new_config = Config()
        new_config.config_data = copy.deepcopy(self.config_data)
        new_config.config_path = self.config_path
        new_config.initialized = True
        return new_config

    def local_clone(self, updates=None):
        """
        Create a new configuration object with updated values.
        :param updates:
        :return:
        """
        if not self.initialized:
            raise RuntimeError("Configuration not loaded. Call load_config first.")
        new_config = self.detach()
        if updates:
            for key, value in updates.items():
                new_config.set(key, value)
        return new_config

    def __getitem__(self, key):
        if not self.initialized:
            raise RuntimeError("Configuration not loaded. Call load_config first.")
        return self.config_data[key]

    def __contains__(self, key):
        if not self.initialized:
            raise RuntimeError("Configuration not loaded. Call load_config first.")
        return key in self.config_data

    def __str__(self):
        status = "Initialized" if self.initialized else "Not Initialized"
        return f"Config Object ({status}):\n  Config Path: {self.config_path}\n  Config Data: {self.config_data}"

    def __repr__(self):
        return f"Config(initialized={self.initialized}, config_path='{self.config_path}')"

