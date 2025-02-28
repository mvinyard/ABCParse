#!/usr/bin/env python

# -- import packages: ---------------------------------------------------------
import ABCParse
import io
import logging
import os
import tempfile
import unittest


# -- Test class: --------------------------------------------------------------
class TestLogging(unittest.TestCase):
    def setUp(self):
        # Capture output for testing
        self.captured_output = io.StringIO()
        
        # Save original handlers and create a test handler
        self.test_handler = logging.StreamHandler(self.captured_output)
        self.test_handler.setFormatter(logging.Formatter("%(message)s"))
        
        # Store original loggers to restore later
        self.original_loggers = {}
        for name, logger in logging.root.manager.loggerDict.items():
            if isinstance(logger, logging.Logger):
                self.original_loggers[name] = {
                    'handlers': logger.handlers.copy(),
                    'level': logger.level,
                    'propagate': logger.propagate
                }
                logger.handlers = []
                logger.addHandler(self.test_handler)
                logger.setLevel(logging.DEBUG)
        
        # Patch the ABCLogger.__init__ method to use our test handler
        self.original_init = ABCParse.logging.ABCLogger.__init__
        
        def patched_init(self_logger, *args, **kwargs):
            # Call the original init
            self.original_init(self_logger, *args, **kwargs)
            
            # Replace handlers with our test handler
            self_logger.logger.handlers = []
            self_logger.logger.addHandler(self.test_handler)
            self_logger.logger.setLevel(logging.DEBUG)
        
        # Apply the patch
        ABCParse.logging.ABCLogger.__init__ = patched_init
        
        # Also patch the default logger
        ABCParse.logging._default_logger.logger.handlers = []
        ABCParse.logging._default_logger.logger.addHandler(self.test_handler)
        ABCParse.logging._default_logger.logger.setLevel(logging.DEBUG)

    def tearDown(self):
        # Restore the original init method
        ABCParse.logging.ABCLogger.__init__ = self.original_init
        
        # Restore original loggers
        for name, config in self.original_loggers.items():
            if name in logging.root.manager.loggerDict:
                logger = logging.getLogger(name)
                logger.handlers = config['handlers']
                logger.level = config['level']
                logger.propagate = config['propagate']
        
        self.captured_output.close()

    def test_get_logger(self):
        """Test that get_logger returns a properly configured logger."""
        logger = ABCParse.logging.get_logger(name="test_logger", level="debug")
        self.assertEqual(logger.name, "test_logger")
        self.assertEqual(logger.level, "debug")
        
        # Test logging at different levels
        logger.debug("Debug message")
        logger.info("Info message")
        
        output = self.captured_output.getvalue()
        self.assertIn("Debug message", output)
        self.assertIn("Info message", output)

    def test_global_log_functions(self):
        """Test the global logging functions."""
        # Clear any previous output
        self.captured_output.truncate(0)
        self.captured_output.seek(0)
        
        ABCParse.logging.set_global_log_level("debug")
        
        ABCParse.logging.debug("Global debug")
        ABCParse.logging.info("Global info")
        ABCParse.logging.warning("Global warning")
        
        output = self.captured_output.getvalue()
        self.assertIn("Global debug", output)
        self.assertIn("Global info", output)
        self.assertIn("Global warning", output)

    def test_log_to_file(self) -> None:
        """Test logging to a file."""
        # Create a temporary file
        fd, temp_path = tempfile.mkstemp()
        os.close(fd)  # Close the file descriptor immediately
        
        try:
            # Test with our custom logger
            logger = ABCParse.logging.ABCLogger(name="TestFileLogger", file_path=temp_path)
            
            try:
                # Log a test message
                test_message = "This is a test message for file logging"
                logger.info(test_message)
                
                # Ensure the logger is properly closed to flush buffers
                logger.close()
                
                # Read the file and check if the message is there
                with open(temp_path, 'r') as f:
                    content = f.read()
                    assert test_message in content, f"Expected message not found in log file. Content: {content}"
            finally:
                # Ensure logger is closed even if assertions fail
                logger.close()
        finally:
            # Clean up the temporary file
            try:
                os.unlink(temp_path)
            except OSError:
                pass  # File might already be deleted or inaccessible

    def test_abc_parse_logging(self) -> None:
        """Test that ABCParse class properly integrates with the logging system."""
        # Create a temporary file for logging
        fd, temp_path = tempfile.mkstemp()
        os.close(fd)  # Close the file descriptor immediately
        
        try:
            # Create a parser with logging to file
            parser = ABCParse(
                name="TestParser",
                log_level="info",
                log_file=temp_path
            )
            
            try:
                # Test that the parser logs correctly
                parser._logger.info("Test message from ABCParse")
                
                # Ensure logger is closed to flush buffers
                parser._logger.close()
                
                # Verify the message was logged to file
                with open(temp_path, 'r') as f:
                    content = f.read()
                    self.assertIn("Test message from ABCParse", content)
            finally:
                # Ensure logger is closed even if assertions fail
                if hasattr(parser, '_logger'):
                    parser._logger.close()
        finally:
            # Clean up the temporary file
            try:
                os.unlink(temp_path)
            except OSError:
                pass  # File might already be deleted or inaccessible


if __name__ == "__main__":
    unittest.main() 