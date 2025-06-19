# TODO list

## Extension of Phase Space Layout

- [ ] Update `AbstractCoordiantes` by adding a distribution, cuts and a `weight` function
- [ ] Update `phase_space_factor` by introducing `in_ps_fac` and `out_ps_fac`, where
      `ps_fac(psp) = in_ps_fac(psp)*out_ps_fac(psp)` and per default `in_ps_fac(::AbstractPSP{T}) =
one(T)`. Therefore, every implementation just needs to implement `out_ps_fac`.
- [ ] Update `PhotonSphericalSystem` by allowing `FixedOmega()` and `FixedOmegaPrime()`(default)

## Extension of `PerturbativeQED`

- [ ] Update `PerturbativeQED` to allow a certain `AbstractMedium` on construction:
      `PerturbativeQED(medium::AbstractMedium)`
- [ ] Implement `Plasma <: AbstractMedium` for a medium which is described by a DSF
- [ ] Implement `Vacuum <: AbstractMedium` for standard perturbative QED (default)

## Extension of `DifferentialCrossSection`

- [x] Update `DifferentialCrossSection` by allowing keyword arguments to fix coordinates
      (Q: is this possible in a type stable way? -> yes; just needs to be moved from nb to
      src)
- [ ] write unit tests for new `DifferentialCrossSection`
- [ ] optimize `insert_by_name`

## Proposal based event generation

- [ ] Update `VegasProposal` to work directly on `Vector{Event}`
- [ ] Insert `VegasMap` as abstraction layer between `VegasGrid` and `VegasProposal` to
      allow different dcs (one for training and one for sampling)
- [ ] Fix the API using `train!`, `build_events(VegasMap)`, `build_events(VegasProposal)`

## Combined Event Generation

- [ ] Implement `CombinedEventGenerator` which combines single event generators (performing Bernoulli on the total cross-sections and then standard generation with the
      selected process)

## Medium

- [x] Consider implementing formula from Arista, Brandt 1984
- [ ] Consider using u = ombar/qbar instead of qbar (like in Arista,Brandt 1984)
- [ ] consider implementing energy loss function and dielectric function
- [ ] Check if RPA works correctly
