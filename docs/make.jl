"""
Make documentation with Documenter.jl
"""

using Documenter

include(joinpath(dirname(@__FILE__), "../src/FullEphemerisPropagator.jl"))


makedocs(
    clean = false,
    build = dirname(@__FILE__),
	modules  = [FullEphemerisPropagator],
    format   = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "FullEphemerisPropagator.jl",
    # options
    pages = [
		"Home" => "index.md",
        "Tutorials" => Any[
            "Basics" => "basics.md",
        ],
		"API" => Any[
			"Core" => "api/api_core.md",
			# "Problem Constructor" => "api/api_create_sft_problem.md",
			# "Core Routines" => "api/api_core.md",
			# "Sims-Flanagan Transcription" => "api/api_simsflanagan.md",
			# "Plotting" => "api/api_plot.md",
		],
    ],
	# assets = [
    #     "./assets/logo.png",
    # ],
)