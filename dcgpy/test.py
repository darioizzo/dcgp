from __future__ import absolute_import as _ai

import unittest as _ut


class test_kernel(_ut.TestCase):
    def runTest(self):
        self.test_double()
        self.test_gdual_double()
        self.test_gdual_vdouble()

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

        # Construction
        ex = expression(inputs=1,
                        outputs=1,
                        rows=1,
                        cols=6,
                        levels_back=6,
                        arity=2,
                        kernels=kernel_set(["sum", "mul", "div", "diff"])(),
                        n_eph=0,
                        seed=33)

        ex = expression(inputs=1,
                        outputs=1,
                        rows=1,
                        cols=6,
                        levels_back=6,
                        arity=2,
                        kernels=kernel_set(["sum", "mul", "div", "diff"])(),
                        n_eph=2,
                        seed=33)
        # Ephemeral value attributes tests
        self.assertEqual(ex.eph_val, [1, 2])
        self.assertEqual(ex.eph_symb, ["c1", "c2"])
        ex.eph_val = [-0.2, 0.3]
        self.assertEqual(ex.eph_val, [-0.2, 0.3])
        ex.eph_symb = ["d1", "d2"]
        self.assertEqual(ex.eph_symb, ["d1", "d2"])

    def test_gdual_double(self):
        from dcgpy import expression_gdual_double as expression
        from dcgpy import kernel_set_gdual_double as kernel_set
        from pyaudi import gdual_double as gdual

        expression(inputs=1,
                   outputs=1,
                   rows=1,
                   cols=6,
                   levels_back=6,
                   arity=2,
                   kernels=kernel_set(["sum", "mul", "div", "diff"])(),
                   n_eph=0,
                   seed=20)

    def test_gdual_vdouble(self):
        from dcgpy import expression_gdual_vdouble as expression
        from dcgpy import kernel_set_gdual_vdouble as kernel_set
        from pyaudi import gdual_vdouble as gdual

        expression(inputs=1,
                   outputs=1,
                   rows=1,
                   cols=6,
                   levels_back=6,
                   arity=2,
                   kernels=kernel_set(["sum", "mul", "div", "diff"])(),
                   n_eph=0,
                   seed=20)

    def test_loss_double(self):
        from dcgpy import expression_double as expression
        from dcgpy import kernel_set_double as kernel_set
        import numpy as np

        ex = expression(inputs=1,
                        outputs=1,
                        rows=1,
                        cols=6,
                        levels_back=6,
                        arity=2,
                        kernels=kernel_set(["sum", "mul", "div", "diff"])(),
                        n_eph=0,
                        seed=33)
        x = 1.
        loss_list = ex.loss([[x]], [ex([x])], "MSE")
        loss_array = ex.loss(np.array([[x]]), np.array([ex([x])]), "MSE")
        self.assertEqual(loss_list, loss_array)

    def test_loss_gdual_double(self):
        from dcgpy import expression_gdual_double as expression
        from dcgpy import kernel_set_gdual_double as kernel_set
        from pyaudi import gdual_double as gdual
        import numpy as np

        ex = expression(inputs=1,
                        outputs=1,
                        rows=1,
                        cols=6,
                        levels_back=6,
                        arity=2,
                        kernels=kernel_set(["sum", "mul", "div", "diff"])(),
                        n_eph=0,
                        seed=33)
        x = gdual(1., "x", 3)
        loss_list = ex.loss([[x]], [ex([x])], "MSE")
        loss_array = ex.loss(np.array([[x]]), np.array([ex([x])]), "MSE")
        self.assertEqual(loss_list, loss_array)

    def test_loss_gdual_vdouble(self):
        from dcgpy import expression_gdual_vdouble as expression
        from dcgpy import kernel_set_gdual_vdouble as kernel_set
        from pyaudi import gdual_vdouble as gdual
        import numpy as np

        ex = expression(inputs=1,
                        outputs=1,
                        rows=1,
                        cols=6,
                        levels_back=6,
                        arity=2,
                        kernels=kernel_set(["sum", "mul", "div", "diff"])(),
                        n_eph=0,
                        seed=33)
        x = gdual([1., 2.], "x", 3)
        loss_list = ex.loss([[x]], [ex([x])], "MSE")
        loss_array = ex.loss(np.array([[x]]), np.array([ex([x])]), "MSE")
        self.assertEqual(loss_list, loss_array)


class test_symbolic_regression(_ut.TestCase):
    def runTest(self):
        from dcgpy import symbolic_regression, generate_koza_quintic, kernel_set_double
        import pygmo as pg
        X, Y = generate_koza_quintic()
        udp = symbolic_regression(
            points=X,
            labels=Y,
            rows=1,
            cols=20,
            levels_back=21,
            arity=2,
            kernels=kernel_set_double(["sum", "diff"])(),
            n_eph=2,
            multi_objective=False,
            parallel_batches=0)
        prob = pg.problem(udp)
        pop = pg.population(prob, 10)
        udp.pretty(pop.champion_x)
        udp.prettier(pop.champion_x)
        # Unconstrained
        self.assertEqual(prob.get_nc(), 0)
        self.assertEqual(prob.get_nic(), 0)
        # Single objective
        self.assertEqual(prob.get_nobj(), 1)
        # Dimensions
        self.assertEqual(prob.get_nix(), 20 * (2 + 1) + 1)
        self.assertEqual(prob.get_nx(), 2 + prob.get_nix())
        # Has gradient and hessians
        self.assertEqual(prob.has_gradient(), True)
        self.assertEqual(prob.has_hessians(), True)


def run_test_suite():
    """Run the full test suite.
    This function will raise an exception if at least one test fails.
    """
    retval = 0
    suite = _ut.TestLoader().loadTestsFromTestCase(test_kernel)
    suite.addTest(test_kernel_set())
    suite.addTest(test_expression())
    suite.addTest(test_symbolic_regression())

    test_result = _ut.TextTestRunner(verbosity=2).run(suite)
    if len(test_result.failures) > 0 or len(test_result.errors) > 0:
        retval = 1
    if retval != 0:
        raise RuntimeError('One or more tests failed.')
