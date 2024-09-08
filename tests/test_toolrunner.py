import pytest
from unittest.mock import patch, MagicMock
from io import StringIO

from bactria.config import Config
from bactria.tool_runner import ToolRunner

@pytest.fixture
def mock_config():
    config = MagicMock(spec=Config)
    config.get.side_effect = lambda key, default=None: {
        'tool_name': 'mock_tool',
        'log_level': 'INFO',
        'tool_parameters': {}
    }.get(key, default)
    return config

@pytest.fixture
def tool_runner(mock_config):
    return ToolRunner(mock_config)

def test_init(tool_runner):
    assert tool_runner.tool_name == 'mock_tool'
    assert isinstance(tool_runner.parameters, dict)

def test_set_get_parameter(tool_runner):
    tool_runner.set_parameter('key', 'value')
    assert tool_runner.get_parameter('key') == 'value'
    assert tool_runner.get_parameter('nonexistent') is None
    assert tool_runner.get_parameter('nonexistent', 'default') == 'default'

def test_build_command(tool_runner):
    tool_runner.set_parameter('param1', 'value1')
    tool_runner.set_parameter('p', 'v')
    command = tool_runner.build_command()
    assert command == ['mock_tool', '--param1', 'value1', '-p', 'v']


@pytest.mark.parametrize("stdout,stderr,return_code,expected_logs", [
    ("Output line 1\nOutput line 2\n", "", 0,
     ["Running command: mock_tool", "Output line 1", "Output line 2", "Command completed successfully"]),
    ("", "Error: something went wrong\n", 1,
     ["Running command: mock_tool", "Error: something went wrong", "Command failed with return code 1"])
])
def test_run(tool_runner, stdout, stderr, return_code, expected_logs, caplog):
    with patch('subprocess.Popen') as mock_popen:
        mock_process = MagicMock()
        mock_process.stdout = StringIO(stdout)
        mock_process.stderr = StringIO(stderr)
        mock_process.wait.return_value = return_code
        mock_popen.return_value = mock_process

        actual_return_code = tool_runner.run()

        assert actual_return_code == return_code
        for log in expected_logs:
            assert log in caplog.text


@pytest.mark.parametrize("input_line,expected_level", [
    ("DEBUG: test message", "debug"),
    ("INFO: test message", "info"),
    ("WARNING: test message", "warning"),
    ("ERROR: test message", "error"),
    ("Regular message", "info")
])
def test_determine_log_level(tool_runner, input_line, expected_level):
    assert tool_runner._determine_log_level(input_line) == expected_level


@pytest.mark.parametrize("input_line,stream,expected_level", [
    ("INFO: test message", "stderr", "info"),
    ("ERROR: test message", "stderr", "error"),
    ("Regular message", "stdout", "info")
])
def test_log_output(input_line, stream, expected_level):
    with patch('bactria.tool_runner.get_formatted_logger') as mock_get_logger:
        mock_logger = MagicMock()
        mock_get_logger.return_value = mock_logger

        mock_config = MagicMock(spec=Config)
        mock_config.get.side_effect = lambda key, default=None: 'INFO' if key == 'log_level' else default

        tool_runner = ToolRunner(mock_config)
        tool_runner._log_output(input_line, stream)

        getattr(mock_logger, expected_level).assert_called_with(input_line.strip())
