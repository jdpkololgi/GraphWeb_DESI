# Historical knot-model comparison

The requested matched replay is complete in the simulation repository:
[methods, metrics and plots](../../TNG/Illustris/docs/p12a_historical_knot_comparison_20260928.md).
Four frozen P8 checkpoints (July U/G-PATCH and their longer-trained descendants)
were run on the identical exposed ph006 50,000-galaxy sample.

At 0.2<z<0.3, current posterior knot recall at 70% precision is 45.2%, versus
31.7–32.5% for the historical models. Current posterior has the highest average
precision in all four redshift shells. July G-PATCH's higher original hard
recall (50.8% versus45.7%) comes with lower precision (57.4% versus69.1%).
The higher old headline recall therefore does not justify switching models.

This does not resolve current knot incompleteness or establish real-DESI
coverage. ph006 was used for current model selection and is not a new blind
test. June whole-wedge FlowJAX was located but not evaluated under the
incompatible P8 input contract. The provisional Loa VAC remains unchanged.
