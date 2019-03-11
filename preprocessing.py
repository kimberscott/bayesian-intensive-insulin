import pandas
import datetime
import math
from ast import literal_eval
import numpy as np

def save_sequences(csvFilename, outFilename, b = 0.5, d = 3, segmentDays = 30, segmentStart = '', sequenceLength = 15, maxGapHours = 16):
    '''Process a CSV file with one row per observation into a CSV with one row per sequence.

    Usage: save_sequences(csvFilename, b = 0.5, d = 3, segmentDays = 30, segment = 8, 
        sequenceLength = 10, maxGapHours = 16)

    csvFilename: path to CSV file to read in
    outFilename: path to CSV file to output
    b: Bin size for insulin and food history vectors, in hours
    d: Total duration of insulin and food history vectors. Impact of data longer than d hours 
        ago should not continue to change.
    segmentDays: how long a segment of data to use for inference at all, in days
    segmentStart: segment start date in format 'yyyy-mm-dd'
    sequenceLength: how many glucose measurements to bin into one sequence
    maxGapHours: how long a single delT (time between glucose measurements) to accept in a sequence

    For each GLUCOSE measurement we obtain:

    Humalog vector: "New injection info" - include only since last glucose measurement 
        (including this timepoint). Binned: [0-b, 1-2b, 2-3b, ... (D-b)-D, D+].
    Carbs vector: "New carbs info" - include only since last glucose measurement 
        (including this timepoint), as for Humalog.
    delT: time since last glucose measurement, in hours.
    Current basal rate: approximate as most recent Lantus dose.

    Then we group the data together into sequences.

    Collect rows 1 thru s, 2 thru (s+1), etc. and gather each set into a new row, with 
    labels MeasBG, delT, Basal, NewCarbs, NewInsulin. Sequences where there's longer than a 
    12-hour gap in measurement are truncated.
    '''

    # Read in the data
    dfFull = pandas.read_csv(csvFilename)
    dfFull = dfFull.reindex(index=dfFull.index[::-1]) # order earliest - latest

    # Pick out data from a particular segment to use, and get time in hours

    timestamps = dfFull['CorrectTimestamp']
    dfFull['Time'] = [datetime.datetime.strptime(t, "%m/%d/%y %H:%M") for t in timestamps] # datetime object times
    df = dfFull.copy()
    
    if not segmentStart:
        startDate = min(df['Time'])
    else:
        startDate = datetime.datetime.strptime(segmentStart, "%Y-%m-%d")
    endDate   = startDate + datetime.timedelta(segmentDays)

    df = df[df['Time'].between(startDate, endDate)]

    df['Hour'] = [t.timestamp()/3600 for t in df['Time']]
    df['Hour'] = df['Hour'] - min(df['Hour'])

    print('Using {} data points between {} and {}'.format(df.shape[0], startDate.strftime('%Y-%m-%d'), endDate.strftime('%Y-%m-%d')))

    # For each GLUCOSE measurement we need to obtain:

    # Humalog vector: "New injection info" - include only since last glucose measurement (including this timepoint). Binned: [0-b, 1-2b, 2-3b, ... (D-b)-D, D].
    # Carbs vector: "New carbs info" - include only since last glucose measurement (including this timepoint), as for Humalog.
    # delT: time since last glucose measurement, in hours.
    # Current basal rate: approximate as most recent Lantus dose.

    nBins = math.ceil(d/b)

    newInsulinEntries = [] # List of insulin since last BG measurement, each in form (timeInHours, insulinUnits)
    newCarbEntries = [] # List of carbs since last BG measurement, each in form (timeInHours, carbGrams)
    lastTime = 0

    # Get Lantus doses from before this time period, if available. Otherwise don't have info until first dose.
    prevLantus = dfFull.loc[(pandas.notna(dfFull["Lantus"])) & (dfFull["Time"] <= startDate), "Lantus"]
    if prevLantus.empty:
        lastLantus = 0
    else:
        lastLantus = prevLantus.values[-1]

    singleEntries = []

    for (index, row) in df.iterrows(): # Loop through data entries. For each entry...
        if pandas.notna(row['Lantus']): # Store most recent Lantus dose
            lastLantus = row['Lantus']
        if pandas.notna(row['Glucose']): # If glucose is recorded, we'll make a new row   
            now = row['Hour']
            # Use list of insulin (and current value) to construct Humalog vector
            humalog = [sum([h for (t,h) in newInsulinEntries 
                            if now-t >= iBin * b and now-t < (iBin+1) * b]) 
                       for iBin in range(nBins+1)]
            # Make sure last bin contains ALL insulin from > nBins * b hours ago
            humalog[-1] = sum([h for (t,h) in newInsulinEntries if now-t >= nBins * b]) 
            humalog[0] += row['Humalog'] if pandas.notna(row['Humalog']) else 0
        
            # Use list of carbs (and current value) to construct carb vector
            carbs = [sum([h for (t,h) in newCarbEntries 
                            if now-t >= iBin * b and now-t < (iBin+1) * b]) 
                     for iBin in range(nBins + 1)]
            # Make sure last bin contains ALL carbs from > nBins * b hours ago
            carbs[-1] = sum([h for (t,h) in newCarbEntries if now-t >= nBins * b]) 
            carbs[0] += row['MealCarbs'] if pandas.notna(row['MealCarbs']) else 0
        
            # Get time since last glucose measurement, in hours
            delT = row['Hour'] - lastTime
        
            # Create a row with delT, glucose, insulin, carbs, basal=lastLantus
            singleEntries.append({
                'MeasBG': row['Glucose'],
                'delT': delT,
                'NewInsulin': humalog,
                'NewCarbs': carbs,
                'Basal': lastLantus
            })
        
            # Clear the lists of new insulin & carbs, and set last T to this T
            newInsulinEntries = []
            newCarbEntries = []
            lastTime = now
        else: # Otherwise, add (T, H) and/or (T, C) to list of current insulin and/or carbs.
            if pandas.notna(row['Humalog']):
                newInsulinEntries.append((row['Hour'], row['Humalog']))
            if pandas.notna(row['MealCarbs']):
                newCarbEntries.append((row['Hour'], row['MealCarbs']))
        
    singleEntries = pandas.DataFrame(singleEntries)

    # Then we need to group the data together into sequences.

    # Collect rows 1 thru s, 2 thru (s+1), etc. and gather each set into a new row, with 
    # appropriate labels. Truncate sequences where there's longer than a 12-hour gap in 
    # measurement. Output a spreadsheet.

    sequences = []

    for iSeq in range(len(singleEntries) - sequenceLength):
        s = singleEntries[iSeq:(iSeq+sequenceLength)].reset_index(drop=True)
        dels = s['delT'].tolist()
        basals = s['Basal'].tolist()
        if all(d <= maxGapHours for d in dels) and all(b > 0 for b in basals):
            sequences.append({
                'MeasBG': s['MeasBG'].tolist(),
                'delT': s['delT'].tolist(),
                'Basal': s['Basal'].tolist(),
                'NewCarbs': s['NewCarbs'].tolist(),
                'NewInsulin': s['NewInsulin'].tolist()
            })
    
    sequences = pandas.DataFrame(sequences)
    sequences.to_csv(outFilename, index=False)
  
def get_protocols(protocolCSV, startDate, endDate=''):
    '''Fetch insulin protocol values from spreadsheet for date or date range.
    
    Usage: 
    get_protocol(protocolCSV, date)
    
        protocolCSV: path to protocol CSV, see protocol.csv for example format
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
    
    get_protocol(protocolCSV, startDate, endDate)
    
        startDate: start of interval for which we want to know protocol ranges
        endDate: end of interval for which we want to know protocol ranges
        
        returns dict with keys:
        CarbRatioBreakfast	
        CarbRatioLunch	
        CarbRatioSnack	
        CarbRatioDinner	
        CarbRatioBedtime	
        CorrFactDay	
        CorrFactNight	
        Lantus
        
        each is a tuple (minValue, maxValue) giving the protocol range within this date range.'''
        
    # Read in the protocol data
    df = pandas.read_csv(protocolCSV)
    
    protocolStarts = df['DateStart']
    protocolStartsFormatted = [datetime.datetime.strftime(datetime.datetime.strptime(pStart, '%m/%d/%y'), '%Y-%m-%d') for pStart in protocolStarts]

    # Get indices of all protocols
    returnSingleValue = not endDate
    if returnSingleValue:
        endDate = startDate
    protocolInds = []
    for iStart in range(len(protocolStartsFormatted)):
        if protocolStartsFormatted[iStart] <= endDate and (iStart == len(protocolStartsFormatted) - 1 or protocolStartsFormatted[iStart + 1] > startDate):
            protocolInds.append(iStart)
    relevantProtocols = df.iloc[protocolInds, :]
    return relevantProtocols

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
            

def scaleByTimeBin(amounts, actionCurve):
    nSeq = len(amounts)
    nBins = len(amounts[0]) 
    return [[amounts[j][k] * actionCurve[k] for k in range(nBins)] for j in range(nSeq)]
    
def extract_amounts_on_board(df, binSize):

    for col in df.columns:
        df.loc[:,col] = df.loc[:,col].apply(lambda x: literal_eval(x))

    nSeg = len(df['Basal'][0])
    nSeq = len(df['Basal'])
    nBins = len(df['NewCarbs'][0][0])

    basal = [ [df['Basal'][i][j] for i in range(nSeq)]
             for j in range(nSeg) ] # basal[i] is array of all basals at first measurement in sequence
    smbg  = [ [df['MeasBG'][i][j] for i in range(nSeq)]
             for j in range(nSeg) ]
    newCarbs = [ [df['NewCarbs'][i][j] for i in range(nSeq)]
             for j in range(nSeg) ] # newCarbs[i] is array of all new-carb arrays at first measurement in sequence. 
    newInsulin = [ [df['NewInsulin'][i][j] for i in range(nSeq)]
             for j in range(nSeg) ]
    delT = [ [df['delT'][i][j] for i in range(nSeq)]
             for j in range(nSeg) ]

    insActCurve = np.multiply(range(0, nBins), 1./(nBins-1)) # TODO: realistic curves; allow scaling
    carbActCurve = np.multiply(range(0, nBins), 1./(nBins-1)) # TODO: realistic curves; allow scaling
    timeBinsDown = np.array(list(reversed(np.multiply(range(1, nBins+1), binSize))))

    # Store the insulin on board and carbs on board that have been added at each timepoint, since last measurement
    newIOB = [scaleByTimeBin(newInsulin[i], list(reversed(insActCurve))) for i in range(nSeg)]
    newCOB = [scaleByTimeBin(newCarbs[i], list(reversed(carbActCurve))) for i in range(nSeg)]

    # Compute the current insulin on board, based on how much remaining IOB is present from last measurement
    # and what IOB has been added since last measurement.

    currentIOB = []
    currentIOB.append(np.array(newIOB[0]))

    for i in range(1, nSeg):
         
        # How much IOB from previous IOB is still there?
        remainingPreviousIOB = [currentIOB[i-1][j] * np.divide(np.maximum(timeBinsDown - delT[i][j], 0), timeBinsDown) 
                                for j in range(nSeq)]   # nSeq x nBins
    
        # Shift it to keep track of how long each dose has left to act
        binsToShift = np.around(np.divide(delT[i], binSize)).astype('int') # len nBins
    
        shiftedPreviousIOB = [np.append([0] * min(binsToShift[j], nBins), prev[:-min(binsToShift[j], nBins)])
                              for (j, prev) in enumerate(remainingPreviousIOB)]
    
        iob = newIOB[i] + shiftedPreviousIOB
    
        currentIOB.append(iob)

    # Analogously, compute current carbs on board 
    # TODO: generalize

    currentCOB = []
    currentCOB.append(np.array(newCOB[0]))

    for i in range(1, nSeg):
         
        # How much COB from previous COB is still there?
        remainingPreviousCOB = [currentCOB[i-1][j] * np.divide(np.maximum(timeBinsDown - delT[i][j], 0), timeBinsDown) 
                                for j in range(nSeq)]   # nSeq x nBins
    
        # Shift it to keep track of how long each food has left to act
        binsToShift = np.around(np.divide(delT[i], binSize)).astype('int') # len nBins
    
        shiftedPreviousCOB = [np.append([0] * min(binsToShift[j], nBins), prev[:-min(binsToShift[j], nBins)])
                              #if binsToShift[j] < nBins else [0] * nBins
                              for (j, prev) in enumerate(remainingPreviousCOB)]
    
        cob = newCOB[i] + shiftedPreviousCOB
    
        currentCOB.append(cob)
    
    iobImpact = [np.minimum(1, np.outer(delT[i], 1/timeBinsDown)) for i in range(nSeg)]
    cobImpact = [np.minimum(1, np.outer(delT[i], 1/timeBinsDown)) for i in range(nSeg)]
    
    return (basal, smbg, newCarbs, newInsulin, delT, insActCurve, carbActCurve,
        currentIOB, currentCOB, iobImpact, cobImpact, nSeq, nSeg)