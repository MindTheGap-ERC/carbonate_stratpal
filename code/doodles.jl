    
    using CarboKitten
    using Unitful
    using CairoMakie
    using CarboKitten.Visualization

    fig = summary_plot("data/init_carbonate_stratpal_1.h5")
    save("figs/test.png", fig)