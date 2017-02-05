""" All logging goes through here """
module Log
export debug, info, warn, error, log, configure, set_log_level

using DocStringExtensions
import Lumberjack

crystal_log = Lumberjack.LumberMill()

for name in [:debug, :info, :warn, :error]
    @eval begin
        $name(msg::AbstractString; kwargs...) =
            Lumberjack.$name(crystal_log, msg; kwargs...)
        $name(msg::AbstractString, args::Dict) =
            Lumberjack.$name(crystal_log, msg, args)
        $name(msg::AbstractString...) = Lumberjack.$name(crystal_log, msg...)
    end
end
log(truck::Lumberjack.TimberTruck, l::Dict) = log(crystal_log, truck, l)
log(mode::AbstractString) = log(crystal_log, mode)
log(mode::AbstractString, msg::AbstractString)  = log(crystal_log, mode, msg)
log(mode::Symbol, msg::AbstractString)  = log(crystal_log, string(mode), msg)
log(mode::AbstractString, msg::AbstractString, args::Dict) =
    log(crystal_log, mode, msg)
log(mode::Symbol, msg::AbstractString, args::Dict) =
    log(crystal_log, string(mode), msg)
log(mode::AbstractString, args::Dict) = log(crystal_log, mode, args)
log(mode::Symbol, args::Dict) = log(crystal_log, string(mode), args)

configure(; kwargs...) = Lumberjack.configure(crystal_log; kwargs...)

"""
    $(SIGNATURES)

Modifies log-level of all "trucks" in Crystals logs. The input should be one of "debug",
"info", "warn", "error", from least to most verbose.
"""
function set_log_level(level::AbstractString="error")
    for (name, truck) in crystal_log.timber_trucks
        Lumberjack.configure(truck, mode=level)
    end
end
set_log_level(level::Symbol) = set_log_level(string(level))

# Sets log level to error by default
set_log_level()

end
