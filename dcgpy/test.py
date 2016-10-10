import unittest as _ut

class test_expression_double(_ut.TestCase):

    def test_construction(self):
        from dcgpy import expression_double as expression
        from dcgpy import function_set_double as function_set



def run_test_suite():
    """Run the full test suite.
    This function will raise an exception if at least one test fails.
    """
    retval = 0
    suite_expression_double = _ut.TestLoader().loadTestsFromTestCase(test_expression_double)
    print("\nRunning tests on expression_double")
    test_result = _ut.TextTestRunner(verbosity=2).run(suite_expression_double)
