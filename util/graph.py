#! /usr/bin/env python3

import plotly.graph_objects as go


def save_histogram(source):
	fig = go.Figure()
	for each_pack_of_length in source:
		fig.add_trace(go.Histogram(
			x=each_pack_of_length["data"],
			name=each_pack_of_length["name"],
			opacity=0.5
		))
	fig.update_layout(barmode='overlay')
	fig.update_traces(opacity=0.5)
	fig.update_layout(yaxis_range=[0, 8000])
	fig.update_layout(xaxis_range=[0, 1000])
	fig.update_layout(title="stribution of aligned length of each score setting")
	fig.update_layout(xaxis_title="length of aligned region (exclude Softclipping)")
	fig.update_layout(yaxis_title="Occurence")
	fig.update_layout(legend_title="score mtx")
	fig.update_layout()
	fig.show()

	#fig.write_image("fig1.jpeg", width=2400, height=1800)

def main():
	fig = go.Figure()
	x  = [12, 24, 60, 121]
	y1 = [39, 79, 200, 396]
	y2 = [21.3, 37.3, 69.4, 133.5]
	fig = go.Figure(data=go.Scatter(x=x, y=y2))
	fig.show()



if __name__ == "__main__":
	main()