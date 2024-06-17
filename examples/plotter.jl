using PyCall 

if !isdefined(Main, :plotter) 
    pushfirst!(pyimport("sys")."path", "examples/")
    plotter = pyimport("plotter")
    importlib = pyimport("importlib")
end

importlib.reload(plotter)