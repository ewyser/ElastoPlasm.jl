# The drawing window without Qt: the editing helpers, the `Editor` and the JSON the view reads. Included by
# DrawQMLExt, and on its own by test_draw.jl, so all of it is tested where QML is not installed.
import ElastoPlasm: inpolygon, inside, get_default, get_option, option_label, option_value, free_numeric, multi_select,
                    shown_name, load_case, save_case, check_case, config_dict, BOUNDARY_KINDS, WORKFLOWS, MATERIAL

include("helpers.jl")
include("editor.jl")
include("view.jl")
