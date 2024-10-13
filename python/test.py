import unittest

loader = unittest.TestLoader()
suite = loader.discover('.', '*_testing.py')

if __name__ == '__main__':
    runner = unittest.TextTestRunner(
        verbosity=2
    )
    runner.run(suite)