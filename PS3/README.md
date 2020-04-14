## PS3 Solution
This solution is slightly different from your PS3 (it does not include the metabolite data from Park et al in the bounds) so your estimated fluxes may be a little different.

### How do I do the QA/QC balance check?
To check if the chemical reactions (contained in the ``Reactions.dat`` file) are balanced, issue the command:

  ```jl
    julia > include("CheckBalances.jl")
  ```
This will formulate the Atom matrix ``A``, and then compute the product of ``transpose(A)*S`` where ``S`` denotes the stoichiometric matrix (contained in the file ``Network.net``). Two arrays are produced: 1) we do not consider the boundary species (we have to nothing/from nothing reactions). These reactions will appear unbalanced with all internal reactions balanced; ii) however, when you include the boundary species, all reactions are balanced.

### How do I estimate the fluxes?
To estimate the Urea flux, issue the command:

  ```jl
    julia > include("Solve.jl")
  ```
The ``Solve`` script formulates the constraints into a [Julia Dictionary](https://docs.julialang.org/en/v1/base/collections/#Dictionaries-1) which is passed to the solver code contained in the ``Flux.jl`` file.
The solver returns a bunch of stuff; the ``objective_value`` and ``flux_array`` arguments contain the Urea flux and
the optimal flux distribution, respectively. The optimal flux that I calculated was approximately: 2.2 mmol/gDW-hr.
Note: this solution does not consider the metabolite levels in the bounds, so your solution could be different.   

### Requirements
The ``Solve.jl`` solution script requires the ``GLPK`` package to the FBA problem. See [GLPK](https://github.com/JuliaOpt/GLPK.jl) for details.
