import pandas as pd
from datetime import date
from datetime import timedelta
import sys
from pathlib import Path

def group_cases_by_state(state):
    
    root_dir = Path(__file__).absolute()
    case_file = root_dir.parent / ".." / "case_data/United_States_COVID-19_Cases_and_Deaths_by_State_over_Time.csv"
    case_data = pd.read_csv(case_file)
    

    state_file = root_dir.parent / ".." / "case_data/states.csv"
    state_list = pd.read_csv(state_file)
    state_dict = state_list.to_dict('records')

    current_state = ""

    for record in state_dict:
        if state.lower() == record['State'].lower():
            current_state = record['Code']

    print(current_state)
    is_current_state = case_data['state'] == current_state
    state_data = case_data[is_current_state]
    state_data = state_data.to_dict('records')
    
    for row in state_data:
        temp_date = row['submission_date'].split("/")
        if len(temp_date[1]) < 2:
            temp_date[1] = f"0{temp_date[1]}"
        if len(temp_date[0]) < 2:
            temp_date[0] = f"0{temp_date[0]}"
        row['submission_date'] = f"{temp_date[2]}-{temp_date[0]}-{temp_date[1]}"

    state_data = pd.DataFrame.from_dict(state_data)
    print(state_data)

    start = date.fromisoformat('2020-01-01')
    end = date.fromisoformat(sorted(state_data.submission_date)[-1])
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
        for row in state_data.itertuples():
            case_date = str(row.submission_date)
            if beginning <= date.fromisoformat(case_date) <= ending:
                count += row.new_case
        temp_frame = pd.DataFrame({"period": period, "count": count}, index = [0])
        period_frame = period_frame.append(temp_frame, ignore_index=True)


    period_frame.to_csv(f"results/{state}_case_period_counts.csv", index = False)
    state_data.to_csv(f"results/{state}_all_case_data.csv")
    print("Right on, all done")