import pandas as pd
import matplotlib as mpl
import seaborn as sns
import seaborn.objects as so
import matplotlib.pyplot as plt

## top-level styling with seaborn
sns.set_theme(
    context='paper', font_scale=0.75, 
    style='darkgrid', 
    palette='deep',
)

## fonts
mpl.rcParams['font.family'] = 'sans-serif'
# sans-serif fonts
preferred_fonts = [
    'Helvetica',
    'Arial',
    'Computer Modern Sans Serif',
]
mpl.rcParams['font.sans-serif'] = preferred_fonts
# use true-font for editable text in illustrator
# default is type-3
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

## layout engine default to constrained
mpl.rcParams['figure.constrained_layout.use'] = True

## resolution
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['savefig.bbox'] = 'tight'

def main(paths, samples, out_stats, out_fig, thresholds):
    dfs = [read_lengths(p, s) for p, s in zip(paths, samples)]
    df = pd.concat(dfs, ignore_index=True)

    report_stats(df, out_stats)

    plot_distribution(df, out_fig, thresholds)
    
def read_lengths(path, sample):
    df = pd.read_csv(path)
    df["sample"] = sample
    return df

def report_stats(df, out):
    (
        df.groupby("condition")["fragment_length"]
        .aggregate(
            min='min', 
            max='max',
            median='median',
            mean='mean',
            std='std'
        )
    ).to_csv(out)

def plot_distribution(df, out, thresholds):
    f = plt.figure()

    p = (
        so.Plot(data=df)
        .facet(row="sample")
        .add(
            so.Bars(alpha=0.4), 
            so.Hist(stat="proportion", binwidth=10, binrange=(0, 1150)),
            x="fragment_length", color="sample",
        )
        .scale(
            x=so.Continuous().tick(every=250)
        )
        .theme(mpl.rcParams)
        .layout(engine='constrained')
        .label(x='Fragment length (bp)', y='Proportion')
    )
    p.on(f).plot()

    for ax in f.axes:
        ymin, ymax = ax.get_ylim()
        ax.vlines(x=thresholds, ymin=ymin, ymax=ymax, colors='grey', alpha=0.4)
        for th in thresholds:
            ax.annotate(text=th, xy=(th+5, ymax*0.9), color='grey')

    f.savefig(out)

if __name__ == '__main__':
    main(
        snakemake.input.csvs,
        snakemake.params.samples,
        snakemake.output.stats,
        snakemake.output.fig,
        snakemake.params.thresholds
    )