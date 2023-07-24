import pandas
import numpy as np
import pymc3 as pm
import matplotlib.pyplot as plt
import datetime
from itertools import chain
import pickle


from preprocessing import (
    save_sequences,
    extract_amounts_on_board,
    get_protocols,
)

full_csv_name = "./data/alldata.csv"
protocol_csv_name = "./data/protocol.csv"

bin_size_h = 0.5
segment_days = 15
start_date_str = "2016-11-01"
end_date_str = "2017-05-01"
experiment_name = "15days_first6mo"

pickle_file = "./processed/results" + experiment_name + ".pickle"
comparison_filename = "./processed/human_vs_inferred_" + experiment_name + ".png"
start_dates = []
end_dates = []
while start_date_str < end_date_str:
    start_dates.append(start_date_str)
    start_date_str = datetime.datetime.strftime(
        datetime.datetime.strptime(start_date_str, "%Y-%m-%d")
        + datetime.timedelta(segment_days),
        "%Y-%m-%d",
    )
    end_dates.append(start_date_str)

# Examples
get_protocols(protocol_csv_name, "2016-11-01", "2017-02-01")

segment_csv_name = "./processed/segment_start_" + str(start_date_str) + ".csv"
save_sequences(
    full_csv_name,
    segment_csv_name,
    bin_size_h=bin_size_h,
    action_duration_h=3,
    segment_days=segment_days,
    segment_start_date_str=start_date_str,
    sequence_n_points=15,
    max_gap_h=16,
)
df = pandas.read_csv(segment_csv_name)
(
    basal,
    smbg,
    new_carbs,
    new_insulin,
    del_t,
    ins_act_curve,
    carb_act_curve,
    current_iob,
    current_cob,
    iob_impact,
    cob_impact,
    nSeq,
    nSeg,
) = extract_amounts_on_board(df, bin_size_h)

# %%


all_data = []
for i_segment, start_date_str in enumerate(start_dates):
    # Get data about what humans inferred parameters to be (not used for inference, just comparison)
    segment_protocols = get_protocols(
        protocol_csv_name, start_date_str, end_dates[i_segment]
    )
    human_cfs = list(
        set(
            chain(
                *[
                    segment_protocols[label].tolist()
                    for label in ["CorrFactDay", "CorrFactNight"]
                ]
            )
        )
    )
    human_crs = list(
        set(
            chain(
                *[
                    segment_protocols[label].tolist()
                    for label in [
                        "CarbRatioBreakfast",
                        "CarbRatioLunch",
                        "CarbRatioSnack",
                        "CarbRatioDinner",
                        "CarbRatioBedtime",
                    ]
                ]
            )
        )
    )
    humanBasals = list(set(segment_protocols["Lantus"].tolist()))

    # Get data about what we actually did
    segment_csv_name = "./processed/segment_start_" + str(start_date_str) + ".csv"
    save_sequences(
        full_csv_name,
        segment_csv_name,
        bin_size_h=bin_size_h,
        action_duration_h=3,
        segment_days=segment_days,
        segment_start_date_str=start_date_str,
        sequence_n_points=15,
        max_gap_h=16,
    )
    df = pandas.read_csv(segment_csv_name)

    (
        basal,
        smbg,
        new_carbs,
        new_insulin,
        del_t,
        ins_act_curve,
        carb_act_curve,
        current_iob,
        current_cob,
        iob_impact,
        cob_impact,
        nSeq,
        nSeg,
    ) = extract_amounts_on_board(df, bin_size_h)

    with pm.Model() as model:
        M = pm.Uniform("Basal", lower=-40, upper=40)  # Prior on basal metabolism
        CF = pm.Normal("CF", mu=250, sd=100)  # Prior on correction factor
        CR = pm.Normal("CR", mu=40, sd=25)  # Prior on correction factor
        sigmaBG = pm.Normal(
            "SigmaBG", mu=50, sd=10
        )  # Prior on measurement error, assuming
        # normally distributed. Note that log-normal distribution of glucose measurement
        # would be more realistic, but a lot slower to run.
        # tauBG = 0.095; # SD of log(BG). If 95% obs w/i 20%, then SD ~ 10%. log(1.1) = 0.095.

        glus = []
        glus.append(pm.Normal("glucose0", mu=150, sd=80, observed=smbg[0]))

        # TODO: infer true IOB, COB

        for i in range(1, nSeg):
            print(i)
            # TODO: use appropriate CR based on when carbs were EATEN. Perhaps unintuitively, will probably
            # be most straightforward to create newCarbsBreakfast, currentCOBBreakfast, cobImpactBreakfast, ...
            # - breaking out into each category - rather than just selecting correct carb ratios here, esp. as
            # we'd need to track which category previously-added carbs fell into (expect carb ratio to
            # depend primarily on when carbs are eaten, not when they are exerting influence on BG)

            # Calculate predicted next glucose value, based on...
            glus.append(
                pm.Normal(
                    "glucose" + str(i),
                    mu=glus[i - 1]
                       + -1  # Last true glucose value
                       * np.dot(new_insulin[i - 1], ins_act_curve)
                       * CF
                       + -1  # Newly-added insulin
                       * np.transpose(
                        [
                            np.dot(current_iob[i - 1][j], iob_impact[i][j])
                            for j in range(nSeq)
                        ]
                    )
                       * CF
                       + np.dot(new_carbs[i - 1], carb_act_curve)  # Insulin on board
                       * CF
                       / CR
                       + np.transpose(  # Newly-added carbs
                        [
                            np.dot(current_cob[i - 1][j], cob_impact[i][j])
                            for j in range(nSeq)
                        ]
                    )
                       * CF
                       / CR
                       + (M - basal[i - 1])  # Carbs on board \
                       * del_t[i - 1]
                       / 24
                       * CF,  # Basal/metabolism mismatch
                    sd=sigmaBG,
                    observed=smbg[i],
                )
            )

        # Inference!
        trace = pm.sample(
            4000, tune=1000, progressbar=True
        )  # draw posterior samples using NUTS sampling

    # plt.figure(figsize=(7, 7))
    # pm.traceplot(trace[100:])
    # plt.tight_layout();

    stats = pm.summary(trace)
    print("Segment beginning ", start_date_str)
    print(stats)
    all_data.append(
        {
            "stats": stats,
            "trace": trace,
            "human_cfs": human_cfs,
            "human_crs": human_crs,
            "humanBasals": humanBasals,
            "start_date_str": start_date_str,
            "end_date_str": end_dates[i_segment],
        }
    )

    with open(pickle_file, "wb") as handle:
        pickle.dump(all_data, handle)

    print("backed up pickled data")

# %%


with open(pickle_file, "rb") as handle:
    all_data = pickle.load(handle)

# %%


fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(11, 11))
plt.xticks(rotation=270)
axs[0].set_title("Correction factors")
axs[0].set_ylabel("mg/dL / unit")
axs[1].set_title("Carb ratios")
axs[1].set_ylabel("g / unit")
axs[2].set_title("Basal dose")
axs[2].set_ylabel("units")

for i_segment in range(len(all_data)):
    data = all_data[i_segment]
    axs[0].plot(
        [data["start_date_str"]] * len(data["human_cfs"]), data["human_cfs"], "ro"
    )
    axs[1].plot(
        [data["start_date_str"]] * len(data["human_crs"]), data["human_crs"], "ro"
    )
    hHuman = axs[2].plot(
        [data["start_date_str"]] * len(data["humanBasals"]), data["humanBasals"], "ro"
    )

    inferred = {
        label: data["stats"].at[label, "mean"] for label in ["CF", "CR", "Basal"]
    }
    inferredError = {
        label: [
            inferred[label] - data["stats"].at[label, "hpd_2.5"],
            data["stats"].at[label, "hpd_97.5"] - inferred[label],
        ]
        for label in ["CF", "CR", "Basal"]
    }

    axs[0].errorbar(
        [data["start_date_str"]], [inferred["CF"]], yerr=[inferredError["CF"]], fmt="bo"
    )
    axs[1].errorbar(
        [data["start_date_str"]], [inferred["CR"]], yerr=[inferredError["CR"]], fmt="bo"
    )
    hInferred = axs[2].errorbar(
        [data["start_date_str"]],
        [inferred["Basal"]],
        yerr=[inferredError["Basal"]],
        fmt="bo",
    )

fig.legend((hHuman[0], hInferred), ("Protocol sheet", "Inferred value"), loc=1)

fig.savefig(comparison_filename)
