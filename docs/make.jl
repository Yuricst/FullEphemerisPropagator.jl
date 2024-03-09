"""
Make documentation with Documenter.jl
"""

using Documenter

include(joinpath(@__DIR__, "../src/FullEphemerisPropagator.jl"))


makedocs(
	modules  = [FullEphemerisPropagator],
    format   = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "FullEphemerisPropagator.jl",
    # options
    pages = [
		"Home" => "index.md",
        "Basics" => "basic.md",
		# "Basics" => Any[
		# 	"Basic Example" => "tutorial.md",
		# 	"Optimization Problem Options" => "optimization_problem.md",
		# 	"Export & Import Solution" => "handle_io.md",
		# 	"Thruster Models" => "thruster.md",
		# ],
		# "Advanced Configurations" => Any[
		# 	"Variable number of nodes" => "advanced_config/variable_nsegs.md",
		# 	"Match-point Definitions" => "advanced_config/matchpoint.md",
		# ],
		# "Advanced Examples" => Any[
		# 	"Mission to asteroids" => "asteroid_missions.md",
		# 	"Global Optimization with MBH" => "use_mbh.md",
		# 	"Arrival to Manifold" => "manifold_arrival.md",
		# ],
		"API" => Any[
			"Core" => "api/api_core.md",
			# "Problem Constructor" => "api/api_create_sft_problem.md",
			# "Core Routines" => "api/api_core.md",
			# "Sims-Flanagan Transcription" => "api/api_simsflanagan.md",
			# "Plotting" => "api/api_plot.md",
		],
    ]
	# assets=[
    #     "assets/psyche_symbol.png",
    # ],
)

deploydocs(
    repo   = "github.com/Yuricst/FullEphemerisPropagator",
    target = "build",
    deps   = nothing,
    make   = nothing,
    push_preview = true
)