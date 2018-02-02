include("c2.jl")

# Reactions
addchangeovertask(n,:Rx1,[:Rx3,:Rx4,:Rx5],[2,1,1])
addchangeovertask(n,:Rx2,[:Rx3,:Rx4,:Rx5],[2,1,1])
addchangeovertask(n,:Rx3,[:Rx1,:Rx2,:Rx4,:Rx5],[1,1,2,1])
addchangeovertask(n,:Rx4,[:Rx3,:Rx5],[2,1])
addchangeovertask(n,:Rx5,[:Rx1,:Rx2,:Rx3,:Rx4],[1,1,1,1])

# Separations
addchangeovertask(n,:Sep1,[:Sep3,:Sep4,:Sep5],[2,2,1])
addchangeovertask(n,:Sep2,[:Sep3,:Sep4,:Sep5],[2,2,1])
addchangeovertask(n,:Sep3,[:Sep1,:Sep2,:Sep5],[1,1,1])
addchangeovertask(n,:Sep4,[:Sep1,:Sep2,:Sep5],[1,1,1])
addchangeovertask(n,:Sep5,[:Sep1,:Sep2,:Sep3,:Sep4],[1,1,1,1])
