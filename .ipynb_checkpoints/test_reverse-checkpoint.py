import numpy as np
import pandas
import hvplot.pandas
from holoviews import opts
import holoviews as hv
import matplotlib

hv.extension('matplotlib')


x_values = np.arange(0, 10)
x_axis = "x"
out_par = "y"
y_values = np.random.randint(1, 20, size=30)
df_group = pandas.DataFrame({
    x_axis: np.append(x_values, np.append(x_values, x_values)),
    out_par: y_values})

print(df_group)
min_vals = []
max_vals = []
percentiles = []
percentile = [25, 50, 75]
for val in x_values:
    df_tmp2 = df_group.loc[df_group[x_axis] == val]
    min_vals.append(df_tmp2[out_par].min())
    max_vals.append(df_tmp2[out_par].max())

    percentiles.append(np.percentile(
        df_tmp2[out_par], percentile, axis=0))
print("Percentiles:\n{}".format(percentiles))
pandas_df = {
    x_axis: x_values,
    "Min": min_vals,
    "Max": max_vals}
percentiles = np.transpose(percentiles)
print("Transposed percentiles:\n{}".format(percentiles))
p_list = []
for i, perc in enumerate(percentile):
    pandas_df["{}. Percentile".format(perc)] = percentiles[i]
    p_list.append("{}. Percentile".format(perc))
pandas_df = pandas.DataFrame(pandas_df)
print(pandas_df)

test_vals = [y_values[0], y_values[10], y_values[20]]
print("######TEST")
print(test_vals)
print("25: {}".format(np.percentile(test_vals, 25)))
print("50: {}".format(np.percentile(test_vals, 50)))
print("75: {}".format(np.percentile(test_vals, 75)))

param_plot = pandas_df.hvplot.line( x=x_axis,
                                    y=p_list,
                                    value_label=out_par,
                                    alpha=0.3)
renderer = hv.Store.renderers['matplotlib'].instance(fig='png', dpi=300)
renderer.save(param_plot, "pics/test.png")