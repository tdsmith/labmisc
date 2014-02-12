# Sample frame documentation

The `sample_frame` parameter accepts an instance of `pandas.DataFrame` with a certain column structure. Unexpected columns are ignored. Column names are case-sensitive and have no synonyms.

Mandatory columns are:
 * Target: a string containing the name of the gene being assayed (e.g. "Actb", "Gapdh")
 * Sample: an object (like a string) uniquely identifying your sample
 * Cq: a floating-point cycle threshold value. Your software may call this Ct or something else; you should rename this column to `Cq`. You can do this in Python using syntax like `my_dataframe.rename(columns={'Ct': 'Cq'}, inplace=True)`.

Optional columns are:
 * Content: recognized values are `'Unknown'` and `'NTC'`. If the column is missing, all wells are assumed to be unknowns. A Content column containing other or missing values will produce an error.
