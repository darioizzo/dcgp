.. _changelog:

Changelog
=========

1.5.1 (Unreleased)
------------------

* All UDAs, :class:`dcgpy.es4cgp`, :class:`dcgpy.mes4cgp`, 
  :class:`dcgpy.moes4cgp`, :class:`dcgpy.momes4cgp` use random mutations rather than active mutations
  to evolve expressions. This improves significantly the algorithm behaviour.
* Bug fix when using a large number of ephemeral constants ... now udas crash is avoided
* New :func:`dcgpy.enable_threading` and :func:`dcgpy.disable_threading` utilities 
  allowing to switch off multithreding entirely. 

1.5 (6/7/2020)
-------------------

New 
~~~

* A new UDA is introduced to solve symbolic regression problems. 
  Its called moes (Multi-Objective Evolutionary Startegy) and completes the 
  evolutionary approaches in dcgp which can now be selected to be memetic or not
  and single objective or multi-objective.

Changes
~~~~~~~

- BREAKING: the API has been made uniform for the four UDAs: :class:`dcgpy.es4cgp`, :class:`dcgpy.mes4cgp`, 
  :class:`dcgpy.moes4cgp`, :class:`dcgpy.momes4cgp` as well as the mutation mechanism. 
  Named parameters have thus changed and default values too. Note that, for example, what
  was *n_mut* in some algos, is now *max_mut*.

- The underlying computations of the symbolic regression optimization problem (UDP) 
  is now performed by obake using a vectorized type. Speed improvements are observed
  of magnitudes between x4 and x100.

- The problem on nans appearing and exceptions being thrown has been solved 
  for :class:`dcgpy.symbolic_regression` by guarding against symengine exceptions
  and by discarding zero columns and rows when inverting hessians for the Newton step of memetic algorithms.

- The UDA :class:`dcgpy.es4cgp` is no longer using a thread bfe to compute the loss. This avoids crashes when pythonic, 
  non thread-safe kernels are used. A bfe can still be set by the user (deprecated in python) after
  the UDA has been instantiated.
  
- Documentation has been improved and all tutorials and examples updated to the new API.

