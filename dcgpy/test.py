from __future__ import absolute_import as _ai

import unittest as _ut

class test_kernel(_ut.TestCase):

    def my_sum(self, x):
        return sum(x)

    def print_my_sum(self, x):
        return "(" + "+".join(x) + ")"

    def test_double(self):
        from dcgpy import kernel_double as kernel

        my_kernel = kernel(self.my_sum, self.print_my_sum, "my_sum_kernel")

        self.assertEqual(my_kernel.__repr__(), "my_sum_kernel")
        self.assertEqual(my_kernel([1, 2, 3]), 6)
        self.assertEqual(my_kernel(["x", "y"]), "(x+y)")

    def test_gdual_double(self):
        from dcgpy import kernel_gdual_double as kernel
        from pyaudi import gdual_double as gdual

        my_kernel = kernel(self.my_sum, self.print_my_sum, "my_sum_kernel")

        self.assertEqual(my_kernel.__repr__(), "my_sum_kernel")
        x = gdual(1, "x", 2)
        y = gdual(2, "y", 2)
        z = gdual(3, "z", 2)
        self.assertEqual(my_kernel([x, y, z]), x + y + z)
        self.assertEqual(my_kernel(["x", "y"]), "(x+y)")

    def test_gdual_vdouble(self):
        from dcgpy import kernel_gdual_vdouble as kernel
        from pyaudi import gdual_vdouble as gdual

        my_kernel = kernel(self.my_sum, self.print_my_sum, "my_sum_kernel")

        self.assertEqual(my_kernel.__repr__(), "my_sum_kernel")
        x = gdual([1, -1], "x", 2)
        y = gdual([2, -2], "y", 2)
        z = gdual([-2, 1], "z", 2)
        self.assertEqual(my_kernel([x, y, z]), x + y + z)
        self.assertEqual(my_kernel(["x", "y"]), "(x+y)")


class test_kernel_set(_ut.TestCase):
    def runTest(self):
            self.test_double()
            self.test_gdual_double()
            self.test_gdual_vdouble()

    def my_sum(self, x):
        return sum(x)

    def print_my_sum(self, x):
        return "(" + "+".join(x) + ")"

    def test_double(self):
        from dcgpy import kernel_set_double as kernel_set
        from dcgpy import kernel_double as kernel
        a = kernel_set(["diff"])
        a.push_back("mul")
        my_kernel = kernel(self.my_sum, self.print_my_sum, "my_sum_kernel")
        a.push_back(my_kernel)
        self.assertEqual(a.__repr__(), "[diff, mul, my_sum_kernel]")
        x = 1
        y = 2
        z = 3
        self.assertEqual(a[0]([x, y, z]), x-y-z)
        self.assertEqual(a[1]([x, y, z]), x*y*z)
        self.assertEqual(a[2]([x, y, z]), x+y+z)

    def test_gdual_double(self):
        from dcgpy import kernel_set_gdual_double as kernel_set
        from dcgpy import kernel_gdual_double as kernel
        from pyaudi import gdual_double as gdual

        a = kernel_set(["diff"])
        a.push_back("mul")
        my_kernel = kernel(self.my_sum, self.print_my_sum, "my_sum_kernel")
        a.push_back(my_kernel)
        self.assertEqual(a.__repr__(), "[diff, mul, my_sum_kernel]")
        x = gdual(1, "x", 2)
        y = gdual(2, "y", 2)
        z = gdual(3, "z", 2)
        self.assertEqual(a[0]([x, y, z]), x-y-z)
        self.assertEqual(a[1]([x, y, z]), x*y*z)
        self.assertEqual(a[2]([x, y, z]), x+y+z)

    def test_gdual_vdouble(self):
        from dcgpy import kernel_set_gdual_vdouble as kernel_set
        from dcgpy import kernel_gdual_vdouble as kernel
        from pyaudi import gdual_vdouble as gdual

        a = kernel_set(["diff"])
        a.push_back("mul")
        my_kernel = kernel(self.my_sum, self.print_my_sum, "my_sum_kernel")
        a.push_back(my_kernel)
        self.assertEqual(a.__repr__(), "[diff, mul, my_sum_kernel]")
        x = gdual([1, -1], "x", 2)
        y = gdual([2, -2], "y", 2)
        z = gdual([-2, 1], "z", 2)
        self.assertEqual(a[0]([x, y, z]), x-y-z)
        self.assertEqual(a[1]([x, y, z]), x*y*z)
        self.assertEqual(a[2]([x, y, z]), x+y+z)


class test_expression(_ut.TestCase):

    def runTest(self):
            self.test_double()
            self.test_gdual_double()
            self.test_gdual_vdouble()
            self.test_loss_double()
            self.test_loss_gdual_double()
            self.test_loss_gdual_vdouble()

    def test_double(self):
        from dcgpy import expression_double as expression
        from dcgpy import kernel_set_double as kernel_set

        ex = expression(1, 1, 1, 6, 6, 2, kernel_set(
            ["sum", "mul", "div", "diff"])(), 32)
        self.assertEqual(ex([1.]), [1])
        self.assertEqual(ex([2.]), [1])
        self.assertEqual(ex([-1.]), [1])
        self.assertEqual(ex([-2.]), [1])

    def test_gdual_double(self):
        from dcgpy import expression_gdual_double as expression
        from dcgpy import kernel_set_gdual_double as kernel_set
        from pyaudi import gdual_double as gdual

        ex = expression(1, 1, 1, 6, 6, 2, kernel_set(
            ["sum", "mul", "div", "diff"])(), 32)
        self.assertEqual(ex([gdual(1, "x", 2)]), [gdual(1)])
        self.assertEqual(ex([gdual(2, "x", 2)]), [gdual(1)])
        self.assertEqual(ex([gdual(-1, "x", 2)]), [gdual(1)])
        self.assertEqual(ex([gdual(-2, "x", 2)]), [gdual(1)])

    def test_gdual_vdouble(self):
        from dcgpy import expression_gdual_vdouble as expression
        from dcgpy import kernel_set_gdual_vdouble as kernel_set
        from pyaudi import gdual_vdouble as gdual

        ex = expression(1, 1, 1, 6, 6, 2, kernel_set(
            ["sum", "mul", "div", "diff"])(), 32)
        self.assertEqual(ex([gdual([1, 2, -1, 2], "x", 2)]),
                         [gdual([1, 1, 1, 1])])

    def test_loss_double(self):
        from dcgpy import expression_double as expression
        from dcgpy import kernel_set_double as kernel_set
        import numpy as np

        ex = expression(1, 1, 1, 6, 6, 2, kernel_set(
            ["sum", "mul", "div", "diff"])(), 32)
        x = 1.
        loss_list = ex.loss([[x]], [ex([x])], "MSE")
        loss_array = ex.loss(np.array([[x]]), np.array([ex([x])]), "MSE")
        self.assertEqual(loss_list, loss_array)

    def test_loss_gdual_double(self):
        from dcgpy import expression_gdual_double as expression
        from dcgpy import kernel_set_gdual_double as kernel_set
        from pyaudi import gdual_double as gdual
        import numpy as np

        ex = expression(1, 1, 1, 6, 6, 2, kernel_set(
            ["sum", "mul", "div", "diff"])(), 32)
        x = gdual(1., "x", 3)
        loss_list = ex.loss([[x]], [ex([x])], "MSE")
        loss_array = ex.loss(np.array([[x]]), np.array([ex([x])]), "MSE")
        self.assertEqual(loss_list, loss_array)


    def test_loss_gdual_vdouble(self):
        from dcgpy import expression_gdual_vdouble as expression
        from dcgpy import kernel_set_gdual_vdouble as kernel_set
        from pyaudi import gdual_vdouble as gdual
        import numpy as np

        ex = expression(1, 1, 1, 6, 6, 2, kernel_set(
            ["sum", "mul", "div", "diff"])(), 32)
        x = gdual([1., 2.], "x", 3)
        loss_list = ex.loss([[x]], [ex([x])], "MSE")
        loss_array = ex.loss(np.array([[x]]), np.array([ex([x])]), "MSE")
        self.assertEqual(loss_list, loss_array)


def run_test_suite():
    """Run the full test suite.
    This function will raise an exception if at least one test fails.
    """
    retval = 0
    suite = _ut.TestLoader().loadTestsFromTestCase(test_kernel)
    suite.addTest(test_kernel_set())
    suite.addTest(test_expression())

    test_result = _ut.TextTestRunner(verbosity=2).run(suite)
    if len(test_result.failures) > 0 or len(test_result.errors) > 0:
        retval = 1
    if retval != 0:
        raise RuntimeError('One or more tests failed.')
