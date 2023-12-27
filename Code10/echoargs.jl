#! /usr/local/bin/julia
#
using ArgParse

const ECHO_PORT = 3000
const ECHO_HOST = "localhost"

function parse_commandline()
  s = ArgParseSettings()
  @add_arg_table s begin
  "--server", "-s"
    help = "hostname of the echo server"
    default = ECHO_HOST
  "--port", "-p"
    help = "port number running the service"
    arg_type = Int
    default = ECHO_PORT
  "--quiet", "-q"
    help = "run quietly, i.e. with no server output"
    action = :store_true
  end
 return parse_args(s)
end

pa = parse_commandline()
host = pa["server"]
eport = pa["port"]
vflag = !pa["quiet"]

println("Host: ",host)
println("Port: ",eport)
println("Flag: ",vflag)

