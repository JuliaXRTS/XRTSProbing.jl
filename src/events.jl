
struct Event{PSP<:AbstractPhaseSpacePoint,T<:Number}
    psp::PSP
    weight::T
end

QEDcore.process(e::Event) = process(e.psp)
QEDcore.model(e::Event) = model(e.psp)
QEDcore.phase_space_definition(e::Event) = phase_space_definition(e.psp)

phase_space_point(e::Event) = e.psp
weight(e::Event) = e.weight
