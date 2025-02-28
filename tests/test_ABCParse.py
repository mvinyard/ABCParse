# -- import packages: ----------------------------------------------------------
import ABCParse
import os
import tempfile
import unittest

# -- test classes: -------------------------------------------------------------
class TestABCParse(unittest.TestCase):

    def setUp(self) -> None:
        """Set up a fresh ABCParse instance for each test."""
        self.parser = ABCParse.ABCParse()
        
    def test_init(self) -> None:
        """Test that ABCParse initializes with expected default values."""
        self.assertFalse(self.parser._BUILT)
        self.assertIsInstance(self.parser._logger, ABCParse.logging.ABCLogger)
        
    def test_build(self) -> None:
        """Test the __build__ method initializes internal structures."""
        self.parser.__build__()
        self.assertTrue(self.parser._BUILT)
        self.assertEqual(self.parser._PARAMS, {})
        self.assertIn("self", self.parser._IGNORE)
        self.assertIn("__class__", self.parser._IGNORE)
        self.assertEqual(self.parser._stored_private, [])
        self.assertEqual(self.parser._stored_public, [])
        
    def test_set(self) -> None:
        """Test the __set__ method for setting attributes."""
        self.parser.__build__()
        
        # Test setting a public attribute
        self.parser.__set__("test_key", "test_value", public=["test_key"], private=[])
        self.assertEqual(self.parser._PARAMS["test_key"], "test_value")
        self.assertEqual(self.parser.test_key, "test_value")
        self.assertIn("test_key", self.parser._stored_public)
        
        # Test setting a private attribute
        self.parser.__set__("private_key", "private_value", public=[], private=["private_key"])
        self.assertEqual(self.parser._PARAMS["private_key"], "private_value")
        self.assertEqual(self.parser._private_key, "private_value")
        self.assertIn("private_key", self.parser._stored_private)
        
    def test_set_existing(self) -> None:
        """Test the __set_existing__ method for updating existing attributes."""
        self.parser.__build__()
        
        # Set up initial attributes
        self.parser.__set__("test_key", "test_value", public=["test_key"], private=[])
        self.parser.__set__("private_key", "private_value", public=[], private=["private_key"])
        
        # Update public attribute
        self.parser.__set_existing__("test_key", "updated_value")
        self.assertEqual(self.parser._PARAMS["test_key"], "updated_value")
        self.assertEqual(self.parser.test_key, "updated_value")
        
        # Update private attribute
        self.parser.__set_existing__("private_key", "updated_private")
        self.assertEqual(self.parser._PARAMS["private_key"], "updated_private")
        self.assertEqual(self.parser._private_key, "updated_private")
        
        # Test kwargs update
        self.parser.kwargs = {"initial": "value"}
        self.parser._PARAMS["kwargs"] = {"initial": "value"}
        self.parser.__set_existing__("kwargs", {"new": "kwarg"})
        self.assertEqual(self.parser.kwargs, {"initial": "value", "new": "kwarg"})
        
        # Test args update
        self.parser.args = (1, 2)
        self.parser._PARAMS["args"] = (1, 2)
        self.parser.__set_existing__("args", (3, 4))
        self.assertEqual(self.parser.args, (1, 2, 3, 4))
        
    def test_stored_property(self) -> None:
        """Test the _STORED property returns combined private and public attributes."""
        self.parser.__build__()
        self.parser.__set__("public_key", "public_value", public=["public_key"], private=[])
        self.parser.__set__("private_key", "private_value", public=[], private=["private_key"])
        
        stored = self.parser._STORED
        self.assertIn("private_key", stored)
        self.assertIn("public_key", stored)
        self.assertEqual(len(stored), 2)
        
    def test_setup_inputs(self) -> None:
        """Test the __setup_inputs__ method."""
        kwargs = {"a": 1, "b": 2, "c": 3}
        
        # Test with empty public and private lists
        public, private = self.parser.__setup_inputs__(kwargs, [], [], [])
        self.assertEqual(public, [])
        self.assertEqual(private, [])
        self.assertTrue(self.parser._BUILT)
        
        # Test with non-empty public list
        public, private = self.parser.__setup_inputs__(kwargs, ["a"], [], [])
        self.assertEqual(public, ["a"])
        self.assertEqual(private, ["a", "b", "c"])
        
        # Test with ignore list
        self.parser.__setup_inputs__(kwargs, [], [], ["ignore_me"])
        self.assertIn("ignore_me", self.parser._IGNORE)
        
    def test_parse(self) -> None:
        """Test the __parse__ method for parsing kwargs."""
        # Test basic parsing
        kwargs = {"a": 1, "b": 2, "self": "ignore_me", "__class__": "ignore_me_too"}
        self.parser.__parse__(kwargs)
        
        self.assertEqual(self.parser._PARAMS, {"a": 1, "b": 2})
        self.assertEqual(getattr(self.parser, "_a", None), 1)
        self.assertEqual(getattr(self.parser, "_b", None), 2)
        self.assertFalse(hasattr(self.parser, "self"))
        
        # Test with public and private lists
        parser2 = ABCParse.ABCParse()
        kwargs = {"a": 1, "b": 2, "c": 3}
        parser2.__parse__(kwargs, public=["a"], private=["b"])
        
        self.assertEqual(getattr(parser2, "a", None), 1)
        self.assertEqual(getattr(parser2, "_b", None), 2)
        self.assertEqual(getattr(parser2, "_c", None), 3)  # Default is private if not specified
        
        # Test with ignore list
        parser3 = ABCParse.ABCParse()
        kwargs = {"a": 1, "b": 2, "ignore_me": 3}
        parser3.__parse__(kwargs, ignore=["ignore_me"])
        
        self.assertEqual(getattr(parser3, "_a", None), 1)
        self.assertEqual(getattr(parser3, "_b", None), 2)
        self.assertFalse(hasattr(parser3, "ignore_me"))
        
    def test_update(self) -> None:
        """Test the __update__ method for updating attributes."""
        # Set up initial state
        self.parser.__parse__({"a": 1, "b": 2})
        
        # Update existing and add new
        self.parser.__update__({"a": 10, "c": 3})
        
        self.assertEqual(self.parser._PARAMS, {"a": 10, "b": 2, "c": 3})
        self.assertEqual(getattr(self.parser, "_a", None), 10)
        self.assertEqual(getattr(self.parser, "_b", None), 2)
        self.assertEqual(getattr(self.parser, "_c", None), 3)
        
        # Test with None values (should be ignored)
        self.parser.__update__({"a": None, "d": 4})
        
        self.assertEqual(self.parser._PARAMS, {"a": 10, "b": 2, "c": 3, "d": 4})
        self.assertEqual(getattr(self.parser, "_a", None), 10)  # Should not be updated
        self.assertEqual(getattr(self.parser, "_d", None), 4)
        
        # Test with public and private lists
        self.parser.__update__({"e": 5, "f": 6}, public=["e"], private=["f"])
        
        self.assertEqual(getattr(self.parser, "e", None), 5)
        self.assertEqual(getattr(self.parser, "_f", None), 6)


class TestFunctionKwargs(unittest.TestCase):
    
    def test_function_kwargs_basic(self) -> None:
        """Test basic functionality of function_kwargs."""
        def test_func(a, b, c=None) -> None:
            pass
        
        kwargs = {"a": 1, "b": 2, "d": 3}
        result = ABCParse.function_kwargs(test_func, kwargs)
        
        self.assertEqual(result, {"a": 1, "b": 2})
        self.assertNotIn("d", result)
        
    def test_function_kwargs_with_obj(self) -> None:
        """Test function_kwargs with an object."""
        def test_func(x, y) -> None:
            pass
        
        class TestObj:
            def __init__(self) -> None:
                self.x = 10
                self.y = 20
                self.z = 30
                
        obj = TestObj()
        result = ABCParse.function_kwargs(test_func, obj=obj)
        
        self.assertEqual(result, {"x": 10, "y": 20})
        self.assertNotIn("z", result)
        
    def test_function_kwargs_with_ignore(self) -> None:
        """Test function_kwargs with ignore list."""
        def test_func(a, b, c) -> None:
            pass
        
        kwargs = {"a": 1, "b": 2, "c": 3}
        result = ABCParse.function_kwargs(test_func, kwargs, ignore=["a", "self", "kwargs"])
        
        self.assertEqual(result, {"b": 2, "c": 3})
        self.assertNotIn("a", result)


class TestAsList(unittest.TestCase):
    
    def test_as_list_single_value(self) -> None:
        """Test as_list with a single value."""
        result = ABCParse.as_list(5)
        self.assertEqual(result, [5])
        
    def test_as_list_list_value(self) -> None:
        """Test as_list with a list value."""
        result = ABCParse.as_list([1, 2, 3])
        self.assertEqual(result, [1, 2, 3])
        
    def test_as_list_with_target_type(self) -> None:
        """Test as_list with target type checking."""
        # Should pass
        result = ABCParse.as_list([1, 2, 3], target_type=int)
        self.assertEqual(result, [1, 2, 3])
        
        # Should pass with single value
        result = ABCParse.as_list(1, target_type=int)
        self.assertEqual(result, [1])
        
        # Should pass with multiple target types
        result = ABCParse.as_list([1, "string"], target_type=[int, str])
        self.assertEqual(result, [1, "string"])
        
        # Should fail
        with self.assertRaises(AssertionError):
            ABCParse.as_list([1, "string"], target_type=int)


class TestLogging(unittest.TestCase):
    
    def test_get_logger(self) -> None:
        """Test get_logger creates a logger with the right name."""
        logger = ABCParse.logging.get_logger("test_logger")
        self.assertEqual(logger.name, "test_logger")
        
    def test_logger_level(self) -> None:
        """Test setting logger level."""
        logger = ABCParse.logging.get_logger("test_level_logger", level="debug")
        self.assertEqual(logger.level, "debug")
        
        logger.set_level("warning")
        self.assertEqual(logger.level, "warning")
        
    def test_file_logging(self) -> None:
        """Test logging to a file."""
        with tempfile.NamedTemporaryFile(delete=False) as temp:
            temp_path = temp.name
            
        try:
            logger = ABCParse.logging.get_logger("file_logger", file_path=temp_path)
            logger.info("Test message")
            logger.close()
            
            with open(temp_path, 'r') as f:
                content = f.read()
                self.assertIn("Test message", content)
        finally:
            # Clean up
            if os.path.exists(temp_path):
                os.remove(temp_path)


if __name__ == "__main__":
    unittest.main()
