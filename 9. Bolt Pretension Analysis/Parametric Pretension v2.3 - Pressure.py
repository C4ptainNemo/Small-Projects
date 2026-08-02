import os
import re
import math
from datetime import datetime
# Unsure if ANSYS has the Python packages, or it needs the user to have them installed.

# This script is for automatically setting the piston pressure and bolt pretension, solving the model, and exporting the solutions tabular data as a csv.
# If there needs to be a change to the models (excluding the forces), delete all but the first, make the changes, then duplicate the required number of times.

# Note: Nothing will print to the Shell until after the script is done, since running it locks the GUI.

def main():
    # 1.1 Define list of design points as (Piston Pressure, M6 Bolts, M4 Bolts) in MPa and Newtons.
    # Ensure the number of models matches the number of design points. Will produce an error and end the script if not.
    # The piston pressure will linearly ramp from 0 to the specified value.
    design_points = [(9, 0, 0), (9, 1000, 500), (9, 3000, 1500), (9, 5000, 2500), (9, 7000, 3500), (9, 10000, 5000)]
    #design_points = [(9, 0, 0)]

    # 1.2. Solve model flag. Set to True to solve the model, False to skip solving. Useful for testing.
    solve_model = False

    # 1.3 Define number of steps. Minimum of 2 steps. 
    # The first is to apply pretension to the bolts, the piston for goes from zero at step 2 to the given value at the last step. 
    # This defines the granularity of the results, but also requires more time to solve. However doubling the steps doesn't give 
    # a doubling in solve time, so this is more efficient than individual solutions at each force.
    # If spliting a max force into increments, add 1 for the the number of steps since the first step is for applying preload to the bolts.
    num_time_steps = 19

    # 1.4. Define the column order
    # The script will look for results starting with these prefixes and write them in this order to the csv.
    # It works by checking if the name of each result in the solution starts with the prefix (Capital Letter)(Number)(Full Stop) using Regular Expressions.
    # If so, it extracts the data and and stores it in a dict keyed by prefix. 
    # Then when writing the csv, it iterates through this list of prefixes to write the columns in the given order.
    # Note: If a prefix is not found in the results, that column will be skipped with a warning printed to the console.
    COLUMN_ORDER = [
        "A1.", "A2.", "A3.", "A4.", "A5.",
        "B1.", "B2.", "B3.", "B4.", "B5.", "B6.", "B7.", "B8.",  "B9.",
        "C1.", "C2.", "C3.", "C4.", "C5.", "C6.",
        "D1.", "D2.", "D3.", "D4."]

    # 1.5. Define output csv directory.
    output_directory = r"C:\Users\Adam\Desktop\Custom Calliper Analysis\Custom Calliper - Oversized Bolts\Results" 
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        
    print("\n---------- Starting Script for Pretension Study ----------")
    print("Results saved to: {}\n".format(output_directory))

    # Check if the number of models equals the number of design points.
    if len(design_points) != len(Model.Analyses):
        raise ValueError("Number of Models ({}) does not equal number of design points ({})".format(
            len(Model.Analyses), len(design_points)))

    # 3. Loop through design points, solving, and exporting results for each.
    for i, dp in enumerate(design_points):
        print("----- Starting iteration for Design Point {} -----".format(i))
        
        # Get the models analysis object.
        analysis = Model.Analyses[i]
        
        '''
        # Debug: Print out the name and type of all the children objects of the Model.
        for child in analysis.Children:
            try:
                print("  Name: {}, Type: {}".format(child.Name, type(child)))
                
                # Can set which child to have it's attributes printed too.
                if "Analysis Settings" in child.Name:
                    print("=== {} ({}) ===".format(name, child_type))
                    
                    # Print all attributes
                    for attr in dir(child):
                        if not attr.startswith("_"):  # Skip private/dunder attributes
                            try:
                                val = getattr(child, attr)
                                # Skip methods, only show properties/values
                                if not callable(val):
                                    print("  {}: {}".format(attr, val))
                            except Exception as e:
                                print("  {}: <error: {}>".format(attr, str(e)))
            except:
                print("  Could not read child")
        '''
        
        piston_pressure = dp[0]
        M6_bolt_pretension = dp[1]
        M4_bolt_pretension = dp[2]
        
        Model.Analyses[i].Name = "{}/{}".format(M6_bolt_pretension, M4_bolt_pretension)
        
        print("Piston Pressure: {} MPa, M6 Bolt Pretension: {} N, M4 Bolt Pretension: {} N".format(
            piston_pressure, M6_bolt_pretension, M4_bolt_pretension))
            

        # Set the number of time steps. If unchanged, skip setting it.
        for child in analysis.Children:
            try:
                if "Analysis Settings" in child.Name:
                    if child.NumberOfSteps == num_time_steps:
                        num_steps_changed = False # Flag for if the number of time steps is changed. Used to ensure piston force updates when changing num_time_steps.
                        #print("Number of time steps already set to {}, skipping".format(num_time_steps))
                    else:
                        child.NumberOfSteps = num_time_steps
                        num_steps_changed = True
            except:
                print("Did not set number of time steps")
        
        # Set the pretension and force values in the model for this design point.
        # Since these are defined as a table of values at each step, the time of each step is defined and the value at that time is set.
        # Set the bolt pretensions to be defined by "Lock" for all but the first time step.
        for child in analysis.Children:
            if "M6 Bolt Pretension" in child.Name:
                # If the first value for the pretension tabular data is already the desired value, skip setting it.
                existing_value_str = str(child.Preload.Output.DiscreteValues[0]).split()[0] # Comes as eg. 5000 [N], so just get number and convert to float.
                existing_value_float = float(existing_value_str)
                
                if existing_value_float == M6_bolt_pretension:
                    #print("M6 Bolt Pretension already set to {} N, skipping".format(M6_bolt_pretension))
                    continue

                child.Preload.Inputs[0].DiscreteValues = [Quantity(0, "s")]
                child.Preload.Output.DiscreteValues = [Quantity(M6_bolt_pretension, "N")]
                for i in range(2, num_time_steps + 1):
                    child.SetDefineBy(i, BoltLoadDefineBy.Lock)
                
            elif "M4 Bolt Pretension" in child.Name:
                # If the first value for the pretension tabular data is already the desired value, skip setting it.
                existing_value_str = str(child.Preload.Output.DiscreteValues[0]).split()[0] # Comes as eg. 5000 [N], so just get number and convert to float.
                existing_value_float = float(existing_value_str)
                
                if existing_value_float == M4_bolt_pretension:
                    #print("M4 Bolt Pretension already set to {} N, skipping".format(M4_bolt_pretension))
                    continue

                child.Preload.Inputs[0].DiscreteValues = [Quantity(0, "s")]
                child.Preload.Output.DiscreteValues = [Quantity(M4_bolt_pretension, "N")]
                for i in range(2, num_time_steps + 1):
                    child.SetDefineBy(i, BoltLoadDefineBy.Lock)
                
            # When applying preload, all other forces must be zero while the preload is applied in the first step.
            # Set the first step as zero, then the last step as the intended piston pressure.
            # The second step (i.e. when the force is first applied) is 0.001, since at zero some things don't seem to respond until the 3rd step.
            elif "Piston Pressure" in child.Name:
                # If the last value for the piston force tabular data is already the desired value, skip setting it.
                # However if the number of time steps was set (occurs if it changed), the num_steps_changed falg will force it to be set.
                try:
                    existing_value_str = str(child.Magnitude.Output.DiscreteValues[2]).split()[0] # Comes as eg. 5000 [N], so just get number and convert to float.
                    existing_value_float = float(existing_value_str)
                except:
                    print("Error reading existing pressure value")
                
                if isclose(existing_value_float, piston_pressure, rel_tol=1e-9) and not num_steps_changed:
                    #print("Piston pressure already set to {} MPa, skipping".format(piston_pressure))
                    continue

                # Set the discrete values at the given times, the inbetween times are linearly interpolated.
                child.Magnitude.Inputs[0].DiscreteValues = [Quantity(0, "s"), Quantity(1, "s"), Quantity(num_time_steps, "s")]
                child.Magnitude.Output.DiscreteValues = [Quantity(0, "MPa"), Quantity(0.000001, "MPa"), Quantity(piston_pressure, "MPa")]
                
        
        # 4. Solve the model
        # Check if the model is unsolved, and the solve_model flag is True.
        status = analysis.Solution.Status
        if (status != SolutionStatusType.Done) and solve_model:
            analysis.Solve(True)
        
        # Print the solution status. If not solved don't save a csv.
        status = analysis.Solution.Status
        print("DP{} solution status: {}".format(i, status))
        
        if status != SolutionStatusType.Done:
            print("Error with solution, no csv written")
            continue

        # 5. Write results to a csv
        try:
            solution = analysis.Solution

            '''
            # Debug: print solution children info
            print("Solution type: {}".format(type(solution)))
            print("Solution children:")
            for child in solution.Children:
                try:
                    print("  Name: {}, Type: {}".format(child.Name, type(child)))
                except:
                    print("  Could not read child")
            '''
            
            timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            
            filename = "Result_Summary_PistonForce_{}N_M6Bolt_{}N_M4Bolt_{}N.csv".format(
                piston_pressure, M6_bolt_pretension, M4_bolt_pretension)
            
            csv_path = os.path.join(output_directory, filename)
            
            solution_elapsed_time = solution.ElapsedRunTime
            element_count = analysis.MeshData.ElementCount
            node_count = analysis.MeshData.NodeCount
            
            # Get the column of piston forces "Tabular Data"
            for child in analysis.Children:
                if "Piston Pressure" in child.Name:
                    # Column index in "Tabular Data" based on the result type. Starts at 1 from the row number.
                    tabular_data_column = 4
                    piston_pressure_raw_data = extract_table(child, tabular_data_column)
                    
                    # The cells are structred "0." or "= 400", so use regex to get just the number. Exclude the first 2 and last values.
                    piston_pressure_column = [float(re.sub(r'[^0-9.]', '', x)) for x in piston_pressure_raw_data[2:-1]]
                    #print(piston_pressure_column)
                    
                    break # Only do once, since they're are multiple named the same.
                    
            # Collect results into a dict keyed by prefix
            results = {}  # { "A1.": (name, value, unit), ... }
            
            for child in solution.Children:
                try:
                    name = child.Name
                    child_type = str(child.DataModelObjectCategory) # Used to classify each type of solution
                    
                    # Debug: Prints the result name and type, used to get the string for the child_type.
                    #print(name); print(child_type)
                    
                    # Define the column index for "Tabular Data" based on the result type. Starts at 1 from the row number.
                    # Since different results have differenet columns, this needs to be set manually.
                    if child_type in ("TotalDeformation", "DirectionalDeformation", "EquivalentStress", "BoltPretensionProbe"):
                        tabular_data_column = 4
                    elif child_type in ("ForceReaction"):
                        tabular_data_column = 6
                    elif child_type in ("DeformationProbe"):
                        tabular_data_column = 5
                    elif child_type in ("SolutionInformation", "TreeGroupingFolder"):
                        tabular_data_column = -1 # To throw an error if it tries to use this
                    else:
                        print("Child type {} not supported".format(child_type))
                        continue
                    
                    # Using regex, check if the childs name has the prefix (Capital Letter)(Number)(Full Stop), eg. A1.
                    # If it does, and the prefix is defined in COLUMN_ORDER, save the name, result type, data, and unit in the results dict keyed by the prefix
                    match = re.match(r"([A-Z]\d+\.)", name)
                    if match:
                        prefix = match.group(1)
                        if prefix in COLUMN_ORDER:
                            
                            data = extract_table(child, tabular_data_column)
                                
                            # Exact the result type and unit from the first row of the "Tabular Data".
                            # Since some of the column headers have the result_type in parentheses, use regex to get three parts; before parentheses, in parentheses, in brackets.
                            header = data[0] # Header
                            data = data[1:] # Exclude the header
                            match = re.match(r'^(.+?)\s*(?:\((.+?)\))?\s*\[(.+)\]$', header)
                            full_name, paren_name, unit = match.groups()
                            result_type = paren_name if paren_name else full_name
                            
                            results[prefix] = (name, result_type, data, unit)
                            
                            # Debug 
                            #print(header); print(name); print(result_type); print(unit); print(data); print("\n")
                            
                except Exception as e:
                    print("  Failed on {}: {}".format(child.Name, str(e)))
            
            # 6. Write CSV
            # Each result's 'data' is a list of values (one per load step).
            # Rows: header row (names), result_type row, unit row, then one row per load step.
            with open(csv_path, "w") as f:
                # Determine number of rows
                num_rows = 0
                for prefix in COLUMN_ORDER:
                    if prefix in results:
                        num_rows = len(results[prefix][2])
                        break

                # Build tabular header, result_type, and unit rows
                name_row        = "Piston Pressure"
                result_type_row = "Piston Pressure"
                unit_row        = "MPa"

                ordered_data = []
                for prefix in COLUMN_ORDER:
                    if prefix in results:
                        name, result_type, data, unit = results[prefix]
                        name_row        += ",{}".format(name)
                        result_type_row += ",{}".format(result_type)
                        unit_row        += ",{}".format(unit)
                        ordered_data.append(data)
                    else:
                        print("  Warning: no result found for prefix {}".format(prefix))
                        name_row        += ","
                        result_type_row += ","
                        unit_row        += ","
                        ordered_data.append([""] * num_rows)

                f.write(name_row        + "\n")
                f.write(result_type_row + "\n")
                f.write(unit_row        + "\n")

                # Write one row per load step
                for row_idx in range(num_rows):
                    piston_pressure_value = piston_pressure_column[row_idx] if row_idx < len(piston_pressure_column) else ""
                    value_row = "{}".format(piston_pressure_value)
                    for col_data in ordered_data:
                        value_row += ",{}".format(col_data[row_idx]) if row_idx < len(col_data) else ","
                    f.write(value_row + "\n")
                    
                # Write fixed metadata as key-value pairs
                f.write("\n")
                f.write("DP,{}\n".format(i))
                f.write("Piston Pressure (scalar),{}\n".format(piston_pressure))
                f.write("M6 Pretension,{}\n".format(M6_bolt_pretension))
                f.write("M4 Pretension,{}\n".format(M4_bolt_pretension))
                f.write("Solution Time,{}\n".format(solution_elapsed_time))
                f.write("Mesh Nodes,{}\n".format(node_count))
                f.write("Mesh Elements,{}\n".format(element_count))
                f.write("Solution Date,{}\n".format(timestamp))
                f.write("\n")
        
            print("Saved Result Summary to: {}".format(filename))
            
        except Exception as e:
            print("Failed to export Result Summary: {}".format(str(e)))
            
    print("\n---------- Ending Script for Pretension Study ----------")

# ---------- Helper Functions ----------
def isclose(a, b, rel_tol=1e-9, abs_tol=0.0):
                return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def extract_table(child, column):
    # Get the desired column data from the "Tabular Data" panel in the GUI.
    # This section is based on this article
    # https://simutechgroup.com/resources/blog/using-python-to-export-tabular-data-from-within-ansys-mechanical/
    
    # Activate the child object
    child.Activate()
    # Create list to store each row in
    data = []
    # Get the tabular data panel object
    Pane = ExtAPI.UserInterface.GetPane(MechanicalPanelEnum.TabularData)
    # I don't know what this does
    Con = Pane.ControlUnknown
    
    for row in range(1, Con.RowsCount + 1):
        data.append(Con.cell(row, tabular_data_column).Text)
    
    return data
    
# ---------- Run Main Function ----------
# This allows the helper functions to be in scope even when they are at the end of the script.
main()
    
