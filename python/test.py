import unittest

loader = unittest.TestLoader()
suite = loader.discover('.', '*_test.py')

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite)