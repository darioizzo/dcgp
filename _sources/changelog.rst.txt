.. _changelog:

Changelog
=========

1.5.1 (Unreleased)
------------------

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

- The underlying computations of the symbolic regression optimizationlroblem (UDP) 
  is now performed by obake using a vectorized type. Speed improvements are observed
  of orders nbetween x4 and x100 depending on cases.

- The problem on nans appearing and exceptions being thrown has been solved
  by guarding against symengine exceptions and by discarding zero columns and rows
  when inverting hessians for the Newton step of memetic algorithms.

- BREAKING: the API has been made uniform for the four UDAs: es4cgp, moes4cgp, mes4cgp, momes4cgp
  as well as the mutation mechanism. Named parameters have thus changed and default values too.

- The UDA es4cgp is no longer using a thread bfe to compute the loss. This avoids crashes when pythonic, 
  non thread-safe kernels are used. A bfe can still be set by the user (deprecated in python) after
  the UDA has been instantiated.
  
- Documentation has been improved and all tutorials and examples updated to the new API.

