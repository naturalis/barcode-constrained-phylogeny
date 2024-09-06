import os
import pytest
import yaml
from bactria.config import Config

# Construct the path to the config file
CONFIG_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config', 'config.yaml')


@pytest.fixture
def config():
    """Fixture to create and load a Config object for each test."""
    cfg = Config()
    cfg.load_config(CONFIG_PATH)
    return cfg


def test_load_config():
    """Test that the config file can be loaded successfully."""
    cfg = Config()
    cfg.load_config(CONFIG_PATH)
    assert cfg.initialized
    assert cfg.config_path == CONFIG_PATH


def test_get_existing_key(config):
    """Test getting an existing key from the config."""
    assert config.get('cpu_cores') == 4


def test_get_non_existing_key(config):
    """Test getting a non-existing key from the config."""
    assert config.get('non_existing_key', 'default') == 'default'


def test_set_key(config):
    """Test setting a new key in the config."""
    config.set('new_key', 'new_value')
    assert config.get('new_key') == 'new_value'


def test_detach(config):
    """Test creating a detached copy of the config."""
    detached = config.detach()
    assert detached.config_data == config.config_data
    assert detached is not config


def test_local_clone(config):
    """Test creating a local clone of the config with updates."""
    clone = config.local_clone({'cpu_cores': 8})
    assert clone.get('cpu_cores') == 8
    assert config.get('cpu_cores') == 4


def test_getitem(config):
    """Test accessing config values using dictionary-like syntax."""
    assert config['marker'] == 'COI-5P'


def test_contains(config):
    """Test checking if a key exists in the config."""
    assert 'marker' in config
    assert 'non_existing_key' not in config


def test_str_representation(config):
    """Test the string representation of the config object."""
    str_repr = str(config)
    assert "Config Object (Initialized)" in str_repr
    assert f"Config Path: {CONFIG_PATH}" in str_repr


def test_repr_representation(config):
    """Test the repr representation of the config object."""
    repr_str = repr(config)
    assert f"Config(initialized=True, config_path='{CONFIG_PATH}')" == repr_str


def test_process_relative_paths():
    """Test that relative paths in the config are processed correctly."""
    cfg = Config()
    cfg.load_config(CONFIG_PATH)

    # Check if the paths in file_names are absolute
    for key, value in cfg.get('file_names').items():
        if key.endswith('_file'):
            assert os.path.isabs(value), f"Path for {key} should be absolute"


def test_config_immutability():
    """Test that changes to a cloned config don't affect the original."""
    original = Config()
    original.load_config(CONFIG_PATH)

    clone = original.local_clone({'cpu_cores': 8})

    assert clone.get('cpu_cores') == 8
    assert original.get('cpu_cores') == 4


def test_error_on_uninitialized_config():
    """Test that operations on uninitialized config raise errors."""
    cfg = Config()

    with pytest.raises(RuntimeError):
        cfg.get('some_key')

    with pytest.raises(RuntimeError):
        cfg.set('some_key', 'value')

    with pytest.raises(RuntimeError):
        cfg.detach()

    with pytest.raises(RuntimeError):
        cfg.local_clone()


def test_load_non_existent_config():
    """Test that loading a non-existent config file raises an error."""
    cfg = Config()
    with pytest.raises(FileNotFoundError):
        cfg.load_config('non_existent_config.yaml')


def test_load_invalid_yaml():
    """Test that loading an invalid YAML file raises an error."""
    import tempfile

    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as temp_file:
        temp_file.write("invalid: yaml: file:")

    cfg = Config()
    with pytest.raises(yaml.YAMLError):
        cfg.load_config(temp_file.name)

    os.unlink(temp_file.name)