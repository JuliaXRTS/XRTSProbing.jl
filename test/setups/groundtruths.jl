
_groundtruth_compute(x) = x
_groundtruth_input_validation(x) = (x > 0)
struct TestException <: QEDbase.AbstractInvalidInputException end
function _groundtruth_valid_input_assert(x)
    _groundtruth_input_validation(x) || throw(TestException())
    nothing
end
_transform_to_invalid(x) = -abs(x)
_groundtruth_post_processing(x, y) = x + y
