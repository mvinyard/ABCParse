import unittest
from ..autoparser import AutoParser


class TestAutoParser(unittest.TestCase):

    def test_init(self):
        parser = AutoParser()
        self.assertEqual(parser._PARAMS, {})

    def test_inspect(self):
        def example_function(a, b, c):
            pass

        parser = AutoParser()
        params = parser._inspect(example_function)
        self.assertEqual(params, ['a', 'b', 'c'])

    def test_init_params(self):
        parser = AutoParser()
        self.assertEqual(parser._init_params, ['self'])

    def test_call_params(self):
        parser = AutoParser()
        self.assertEqual(parser._call_params, ['self'])

    def test_parse_params(self):
        parser = AutoParser()
        self.assertEqual(parser._parse_params, ['self', 'kwargs', 'ignore', 'private', 'public', 'kwargs_key'])

    def test_collected_params(self):
        parser = AutoParser()
        expected = ['self', 'self', 'self', 'kwargs', 'ignore', 'private', 'public', 'kwargs_key']
        self.assertEqual(parser._collected_params, expected)

    def test_collect_literal_kwargs(self):
        parser = AutoParser()
        parser._collect_literal_kwargs({'a': 1, 'b': 2})
        self.assertEqual(parser._PARAMS, {'a': 1, 'b': 2})
        self.assertEqual(parser.a, 1)
        self.assertEqual(parser.b, 2)

    def test_hide(self):
        parser = AutoParser()
        self.assertEqual(parser.__hide__('key'), '_key')

    def test_collect(self):
        parser = AutoParser()
        parser.__collect__('c', 3)
        self.assertEqual(parser._PARAMS, {'c': 3})
        self.assertEqual(parser.c, 3)

    def test_parse(self):
        parser = AutoParser()
        locals_dict = {'a': 1, 'b': 2, '_c': 3, 'kwargs': {'d': 4}}
        parser.__parse__(locals_dict, public=['a', 'b'], private=['_c'], kwargs_key='kwargs')
        self.assertEqual(parser._PARAMS, {'a': 1, 'b': 2, '_c': 3, 'd': 4})
        self.assertEqual(parser.a, 1)
        self.assertEqual(parser.b, 2)
        self.assertEqual(parser._c, 3)
        self.assertEqual(parser.d, 4)


if __name__ == "__main__":
    unittest.main()
