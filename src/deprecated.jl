#=
This file gathers all the deprecated names (structures, functions, or variables) which will be removed in the future major release.

- Before the major release, the deprecated names will just throw errors when they are called.
- If the deprecated names were once exported, we will still export them here until next major release.
- If we decide to push a major release, cleanup this file.

Example:

export deprecated_foo

function deprecated_foo(args...; kwargs...)
    error("`deprecated_foo` has been deprecated and will be removed in next major release, please use `new_foo` instead.")
end
=#
