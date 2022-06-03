from itertools import count
import sys
import pandas as pd
from datetime import date
from datetime import timedelta
import sys

case_data = pd.read_csv(sys.argv[1])
case_data = case_data.to_dict('records')

for row in case_data:
    temp_date = row['submission_date'].split("/")
    if len(temp_date[1]) < 2:
        temp_date[1] = f"0{temp_date[1]}"
    if len(temp_date[0]) < 2:
        temp_date[0] = f"0{temp_date[0]}"
    row['submission_date'] = f"{temp_date[2]}-{temp_date[0]}-{temp_date[1]}"

case_data = pd.DataFrame.from_dict(case_data)

start = date.fromisoformat('2020-01-01')
end = date.fromisoformat(sorted(case_data['submission_date'])[-1])
week = timedelta(days=7)

periods = []
while start < end:
    period = f"{start}_{start + week}"
    periods.append(period)
    start = start + week

period_frame = pd.DataFrame(columns = ["period", "count"])
for period in periods:
    beginning = date.fromisoformat(period.split("_")[0])
    ending = date.fromisoformat(period.split("_")[1])
    count = 0
    for row in case_data.itertuples():
        case_date = str(row.submission_date)
        if beginning <= date.fromisoformat(case_date) <= ending:
            count += row.new_case
    temp_frame = pd.DataFrame({"period": period, "count": count}, index = [0])
    period_frame = period_frame.append(temp_frame, ignore_index=True)


period_frame.to_csv(sys.argv[2], index = False)
print("Right on, all done")