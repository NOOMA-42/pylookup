import unittest

if __name__ == "__main__":
    # ellipticcurve
    loader = unittest.TestLoader()
    # util_suite = loader.discover(start_dir='./common_util/algebra', pattern='*_test.py')
    src_suite = loader.discover(start_dir='./src', pattern='*_test.py')
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(src_suite)