

@inline function _select_accepted(mask, trail_sample)
    # use fancy indexing here: only allowed on CPU
    # TODO: is this event based, or weight based?
    # - consider using an Event type and an update mechanism (e.g. using Refs or arrays)
    return trail_sample[mask]
end
