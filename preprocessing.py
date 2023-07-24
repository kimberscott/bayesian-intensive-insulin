import pandas as pd
import datetime
import math
from ast import literal_eval
import numpy as np


def save_sequences(
    csv_filename: str,
    out_filename: str,
    bin_size_h: float = 0.5,
    action_duration_h: float = 3,
    segment_days: float = 30,
    segment_start_date_str: str = "",
    sequence_n_points: int = 15,
    max_gap_h: float = 16,
) -> None:
    """Process a CSV file with one row per observation into a CSV with one row per sequence.

    Usage: save_sequences(csv_filename, bin_size_h = 0.5, d = 3, segment_days = 30, segment = 8,
        sequence_length = 10, maxGapHours = 16)

    :param csv_filename: path to CSV file to read in
    :param out_filename: path to CSV file to output
    :param bin_size_h: Bin size for insulin and food history vectors, in hours
    :param action_duration_h: Total duration of insulin and food history vectors. Impact of data longer than d hours
        ago should not continue to change.
    :param segment_days: how long a segment of data to use for inference at all, in days
    :param segment_start_date_str: segment start date in format 'yyyy-mm-dd'
    :param sequence_n_points: how many glucose measurements to bin into one sequence
    :param max_gap_h: how long a single del_t (time between glucose measurements) to accept in a sequence

    For each GLUCOSE measurement we obtain:

    Humalog vector: "New injection info" - include only since last glucose measurement
        (including this timepoint). Binned: [0-bin_size_h, 1-2b, 2-3b, ... (D-bin_size_h)-D, D+].
    Carbs vector: "New carbs info" - include only since last glucose measurement
        (including this timepoint), as for Humalog.
    del_t: time since last glucose measurement, in hours.
    Current basal rate: approximate as most recent Lantus dose.

    Then we group the data together into sequences.

    Collect rows 1 thru s, 2 thru (s+1), etc. and gather each set into a new row, with
    labels MeasBG, del_t, Basal, NewCarbs, NewInsulin. Sequences where there's longer than a
    12-hour gap in measurement are truncated.
    """

    # Read in the data
    df_full = pd.read_csv(csv_filename)
    df_full = df_full.reindex(index=df_full.index[::-1])  # order earliest - latest

    # Pick out data from a particular segment to use, and get time in hours

    timestamps = df_full["CorrectTimestamp"]
    df_full["Time"] = [
        datetime.datetime.strptime(t, "%m/%d/%y %H:%M") for t in timestamps
    ]  # datetime object times
    df = df_full.copy()

    if not segment_start_date_str:
        start_date = min(df["Time"])
    else:
        start_date = datetime.datetime.strptime(segment_start_date_str, "%Y-%m-%d")
    end_date = start_date + datetime.timedelta(segment_days)

    df = df[df["Time"].between(start_date, end_date)]

    df["Hour"] = [t.timestamp() / 3600 for t in df["Time"]]
    df["Hour"] = df["Hour"] - min(df["Hour"])

    print(
        "Using {} data points between {} and {}".format(
            df.shape[0], start_date.strftime("%Y-%m-%d"), end_date.strftime("%Y-%m-%d")
        )
    )

    # For each GLUCOSE measurement we need to obtain:

    # Humalog vector: "New injection info" - include only since last glucose measurement (including this timepoint).
    #    Binned: [0-bin_size_h, 1-2b, 2-3b, ... (D-bin_size_h)-D, D].
    # Carbs vector: "New carbs info" - include only since last glucose measurement (including this timepoint),
    #    as for Humalog.
    # del_t: time since last glucose measurement, in hours.
    # Current basal rate: approximate as most recent Lantus dose.

    n_bins = math.ceil(action_duration_h / bin_size_h)

    new_insulin_entries = (
        []
    )  # List of insulin since last BG measurement, each in form (timeInHours, insulinUnits)
    new_carb_entries = (
        []
    )  # List of carbs since last BG measurement, each in form (timeInHours, carbGrams)
    last_time = 0

    # Get Lantus doses from before this time period, if available. Otherwise don't have info until first dose.
    prev_lantus = df_full.loc[
        (pd.notna(df_full["Lantus"])) & (df_full["Time"] <= start_date), "Lantus"
    ]
    if prev_lantus.empty:
        last_lantus = 0
    else:
        last_lantus = prev_lantus.values[-1]

    single_entries = []

    for index, row in df.iterrows():  # Loop through data entries. For each entry...
        if pd.notna(row["Lantus"]):  # Store most recent Lantus dose
            last_lantus = row["Lantus"]
        if pd.notna(row["Glucose"]):  # If glucose is recorded, we'll make a new row
            now = row["Hour"]
            # Use list of insulin (and current value) to construct Humalog vector
            humalog = [
                sum(
                    [
                        h
                        for (t, h) in new_insulin_entries
                        if iBin * bin_size_h <= now - t < (iBin + 1) * bin_size_h
                    ]
                )
                for iBin in range(n_bins + 1)
            ]
            # Make sure last bin contains ALL insulin from > n_bins * bin_size_h hours ago
            humalog[-1] = sum(
                [h for (t, h) in new_insulin_entries if now - t >= n_bins * bin_size_h]
            )
            humalog[0] += row["Humalog"] if pd.notna(row["Humalog"]) else 0

            # Use list of carbs (and current value) to construct carb vector
            carbs = [
                sum(
                    [
                        h
                        for (t, h) in new_carb_entries
                        if iBin * bin_size_h <= now - t < (iBin + 1) * bin_size_h
                    ]
                )
                for iBin in range(n_bins + 1)
            ]
            # Make sure last bin contains ALL carbs from > n_bins * bin_size_h hours ago
            carbs[-1] = sum([h for (t, h) in new_carb_entries if now - t >= n_bins * bin_size_h])
            carbs[0] += row["MealCarbs"] if pd.notna(row["MealCarbs"]) else 0

            # Get time since last glucose measurement, in hours
            del_t = row["Hour"] - last_time

            # Create a row with del_t, glucose, insulin, carbs, basal=last_lantus
            single_entries.append(
                {
                    "MeasBG": row["Glucose"],
                    "del_t": del_t,
                    "NewInsulin": humalog,
                    "NewCarbs": carbs,
                    "Basal": last_lantus,
                }
            )

            # Clear the lists of new insulin & carbs, and set last T to this T
            new_insulin_entries = []
            new_carb_entries = []
            last_time = now
        else:  # Otherwise, add (T, H) and/or (T, C) to list of current insulin and/or carbs.
            if pd.notna(row["Humalog"]):
                new_insulin_entries.append((row["Hour"], row["Humalog"]))
            if pd.notna(row["MealCarbs"]):
                new_carb_entries.append((row["Hour"], row["MealCarbs"]))

    single_entries = pd.DataFrame(single_entries)

    # Then we need to group the data together into sequences.

    # Collect rows 1 thru s, 2 thru (s+1), etc. and gather each set into a new row, with
    # appropriate labels. Truncate sequences where there's longer than a 12-hour gap in
    # measurement. Output a spreadsheet.

    sequences = []

    for iSeq in range(len(single_entries) - sequence_n_points):
        s = single_entries[iSeq: (iSeq + sequence_n_points)].reset_index(drop=True)
        dels = s["del_t"].tolist()
        basals = s["Basal"].tolist()
        if all(d <= max_gap_h for d in dels) and all(b > 0 for b in basals):
            sequences.append(
                {
                    "MeasBG": s["MeasBG"].tolist(),
                    "del_t": s["del_t"].tolist(),
                    "Basal": s["Basal"].tolist(),
                    "NewCarbs": s["NewCarbs"].tolist(),
                    "NewInsulin": s["NewInsulin"].tolist(),
                }
            )

    sequences = pd.DataFrame(sequences)
    sequences.to_csv(out_filename, index=False)


def get_protocols(protocol_csv: str, start_date_str: str, end_date_str: str = "") -> pd.DataFrame:
    """Fetch insulin protocol values from spreadsheet for date or date range.

    Usage:
    get_protocol(protocol_csv, date)

        protocol_csv: path to protocol CSV, see protocol.csv for example format
        date: date on which we want to know the protocol, in format 'yyyy-mm-dd'

        returns dict with keys:
        CarbRatioBreakfast
        CarbRatioLunch
        CarbRatioSnack
        CarbRatioDinner
        CarbRatioBedtime
        CorrFactDay
        CorrFactNight
        Lantus

        Each is a single value, giving the protocol on this date.

    get_protocol(protocol_csv, start_date_str, end_date_str)

        start_date_str: start of interval for which we want to know protocol ranges
        end_date_str: end of interval for which we want to know protocol ranges

        returns dict with keys:
        CarbRatioBreakfast
        CarbRatioLunch
        CarbRatioSnack
        CarbRatioDinner
        CarbRatioBedtime
        CorrFactDay
        CorrFactNight
        Lantus

        each is a tuple (minValue, maxValue) giving the protocol range within this date range.
    """

    # Read in the protocol data
    df = pd.read_csv(protocol_csv)

    protocol_starts = df["DateStart"]
    protocol_starts_formatted = [
        datetime.datetime.strftime(
            datetime.datetime.strptime(pStart, "%m/%d/%y"), "%Y-%m-%d"
        )
        for pStart in protocol_starts
    ]

    # Get indices of all protocols
    return_single_value = not end_date_str
    if return_single_value:
        end_date_str = start_date_str
    protocol_inds = []
    for iStart in range(len(protocol_starts_formatted)):
        if protocol_starts_formatted[iStart] <= end_date_str and (
                iStart == len(protocol_starts_formatted) - 1
                or protocol_starts_formatted[iStart + 1] > start_date_str
        ):
            protocol_inds.append(iStart)
    relevant_protocols = df.iloc[protocol_inds, :]
    return relevant_protocols


#     labels = [  "CarbRatioBreakfast",
#                 "CarbRatioLunch",
#                 "CarbRatioSnack",
#                 "CarbRatioDinner",
#                 "CarbRatioBedtime",
#                 "CorrFactDay",
#                 "CorrFactNight",
#                 "Lantus"]
#     if returnSingleValue:
#         assert(len(protocolInds) == 1)
#         return {label: float(relevantProtocols[label].iloc[0]) for label in labels}
#     else:
#         assert(len(protocolInds) >= 1)
#         return {label: (float(min(relevantProtocols[label])),
#                         float(max(relevantProtocols[label]))) for label in labels}


def scale_by_time_bin(amounts, actionCurve):
    n_seq = len(amounts)
    n_bins = len(amounts[0])
    return [[amounts[j][k] * actionCurve[k] for k in range(n_bins)] for j in range(n_seq)]


def extract_amounts_on_board(df, bin_size):
    for col in df.columns:
        df.loc[:, col] = df.loc[:, col].apply(lambda x: literal_eval(x))

    n_seg = len(df["Basal"][0])
    n_seq = len(df["Basal"])
    n_bins = len(df["NewCarbs"][0][0])

    basal = [
        [df["Basal"][i][j] for i in range(n_seq)] for j in range(n_seg)
    ]  # basal[i] is array of all basals at first measurement in sequence
    smbg = [[df["MeasBG"][i][j] for i in range(n_seq)] for j in range(n_seg)]
    new_carbs = [
        [df["NewCarbs"][i][j] for i in range(n_seq)] for j in range(n_seg)
    ]  # new_carbs[i] is array of all new-carb arrays at first measurement in sequence.
    new_insulin = [[df["NewInsulin"][i][j] for i in range(n_seq)] for j in range(n_seg)]
    del_t = [[df["del_t"][i][j] for i in range(n_seq)] for j in range(n_seg)]

    ins_act_curve = np.multiply(
        range(0, n_bins), 1.0 / (n_bins - 1)
    )  # TODO: realistic curves; allow scaling
    carb_act_curve = np.multiply(
        range(0, n_bins), 1.0 / (n_bins - 1)
    )  # TODO: realistic curves; allow scaling
    time_bins_down = np.array(list(reversed(np.multiply(range(1, n_bins + 1), bin_size))))

    # Store the insulin on board and carbs on board that have been added at each timepoint, since last measurement
    new_iob = [
        scale_by_time_bin(new_insulin[i], list(reversed(ins_act_curve)))
        for i in range(n_seg)
    ]
    new_cob = [
        scale_by_time_bin(new_carbs[i], list(reversed(carb_act_curve)))
        for i in range(n_seg)
    ]

    print("new_iob: {} x {} x {}".format(len(new_iob), len(new_iob[0]), len(new_iob[0][0])))

    def get_iob(new_iob):
        """
        Compute the current insulin (/carbs/etc.) on board, based on how much remaining
        IOB is present from last measurement and what IOB has been added since last
        measurement.

        new_iob has dimensions n_seg x n_seq x n_bins
        """

        current_iob = [np.array(new_iob[0])]

        for i in range(1, n_seg):
            # How much IOB from previous IOB is still there?
            remaining_previous_iob = [
                current_iob[i - 1][j]
                * np.divide(np.maximum(time_bins_down - del_t[i][j], 0), time_bins_down)
                for j in range(n_seq)
            ]  # n_seq x n_bins

            # Shift it to keep track of how long each dose has left to act
            bins_to_shift = np.around(np.divide(del_t[i], bin_size)).astype(
                "int"
            )  # len n_bins

            shifted_previous_iob = [
                np.append(
                    [0] * min(bins_to_shift[j], n_bins),
                    prev[: -min(bins_to_shift[j], n_bins)],
                )
                for (j, prev) in enumerate(remaining_previous_iob)
            ]

            iob = new_iob[i] + shifted_previous_iob

            current_iob.append(iob)

        return current_iob

    current_iob = get_iob(new_iob)
    current_cob = get_iob(new_cob)

    iob_impact = [
        np.minimum(1, np.outer(del_t[i], 1 / time_bins_down)) for i in range(n_seg)
    ]
    cob_impact = [
        np.minimum(1, np.outer(del_t[i], 1 / time_bins_down)) for i in range(n_seg)
    ]

    # TODO: fix
    return (
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
        n_seq,
        n_seg,
    )
