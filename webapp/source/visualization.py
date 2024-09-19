import plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd

# def plot_distribution_proportions(sequences, smoothing):
#     fig = go.Figure()
#     for p in sequences:
#         title = p['type']
#         if 'value_counts' in p:
#             df = pd.DataFrame(p['value_counts'])
#             #total_proportion = np.round(np.sum([item['total_proportion'] for item in p['peaks']]),2) if len(p['peaks']) > 0 else 0
#             total_proportion = p['total_proportion']
#             total_reads = p['total_reads']
#             noise = p['noise_level']
#             #
#             # peaks_info = ""
#             # for i,item in enumerate(p['peaks']):
#             #     idx = item['peak_index']
#             #     proportion = item['total_proportion']
#             #     peaks_info += f'{idx}({proportion})  '
#             fig.add_trace(go.Scatter(
#                 x=df['index'],
#                 y=df['proportion'],
#                 mode='lines',
#                 name=f'{title}: Total Proportion = {total_proportion}, number of reads = {total_reads}, level of noise = {noise}',
#                 )
#             )
#             if len(sequences) == 1:
#                 # Add dots at specific x positions
#                 peak_positions = [int(i['peak_index']) for i in p['peaks']]
#                 fig.add_trace(go.Scatter(
#                     x=peak_positions,
#                     y=[df.loc[df['index'] == x, 'proportion'].values[0] for x in peak_positions],
#                     text=[i['total_proportion'] for i in p['peaks']],
#                     textposition='top center',
#                     textfont=dict(color='red'),
#                     mode='markers+text',
#                     marker=dict(size=8, color='red'),
#                     showlegend=False
#                 ))
#                 if smoothing != 'None':
#                     # Smoothing line
#                     x_smooth = df['index']
#                     y_smooth = df['smoothed']
#                     fig.add_trace(go.Scatter(
#                         x=x_smooth,
#                         y=y_smooth,
#                         mode='lines',
#                         showlegend=False,
#                         line=dict(color='red')
#                     ))
#
#     fig.update_layout(
#         width=1200,
#         height=600,
#         legend=dict(
#             x=0,
#             y=1.0,
#             xanchor='left',
#             yanchor='bottom'
#         ),
#         barmode='overlay',
#         title={
#             'text': "Positional distribution of occurrences proportional to the total number of sequences",
#             'y': 0.05,
#             'x': 0.5,
#             'xanchor': 'center',
#             'yanchor': 'bottom'}
#     )
#     fig.update_layout(barmode='overlay')
#     fig.update_traces(opacity=0.75)
#
#     return fig


def plot_distribution_proportions(sequences, smoothing):
    fig = go.Figure()

    for p in sequences:
        title = p['type']
        if 'value_counts' in p:
            df = pd.DataFrame(p['value_counts'])
            total_proportion = p['total_proportion']
            total_reads = p['total_reads']
            noise = p['noise_level']

            # Proportion line on the first Y-axis
            fig.add_trace(go.Scatter(
                x=df['index'],
                y=df['proportion'],
                mode='lines',
                name=f'{title}: Total Proportion = {total_proportion}, number of reads = {total_reads}, level of noise = {noise}',
                yaxis='y1'  # Specify the first Y-axis
            ))

            if len(sequences) == 1:
                # Add dots at specific x positions (for proportions)
                peak_positions = [int(i['peak_index']) for i in p['peaks']]
                fig.add_trace(go.Scatter(
                    x=peak_positions,
                    y=[df.loc[df['index'] == x, 'proportion'].values[0] for x in peak_positions],
                    text=[i['total_proportion'] for i in p['peaks']],
                    textposition='top center',
                    textfont=dict(color='red'),
                    mode='markers+text',
                    marker=dict(size=8, color='red'),
                    showlegend=False
                ))

                # Smoothing line
                if smoothing != 'None':
                    x_smooth = df['index']
                    y_smooth = df['smoothed']
                    fig.add_trace(go.Scatter(
                        x=x_smooth,
                        y=y_smooth,
                        mode='lines',
                        showlegend=False,
                        line=dict(color='red')
                    ))

    # Update layout with dual Y-axis (proportion on left, reads on right)
    fig.update_layout(
        width=1200,
        height=600,
        legend=dict(
            x=0,
            y=1.0,
            xanchor='left',
            yanchor='bottom'
        ),
        barmode='overlay',
        title={
            'text': "Positional distribution of occurrences proportional to the total number of sequences (Proportion) and Reads",
            'y': 0.05,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'bottom'
        },
        # Define the primary Y-axis for Proportion
        yaxis=dict(
            title="Proportion",
            side="left",
            showgrid=True
        ),
        # Define the secondary Y-axis for Reads (without plotting reads)
        yaxis2=dict(
            title="Reads",
            side="right",
            overlaying="y",  # Allows the second Y-axis to overlay on the same plot
            showgrid=False,  # Disable grid for the second axis
            range=[0, max([p['total_reads'] for p in sequences])],  # Set range for the second Y-axis based on reads
            zeroline=True,  # Ensure zero line is shown for visibility
            showline=True,  # Ensure the axis line is visible
            ticks="outside",  # Show ticks outside the plot for the second axis
            tickmode='auto'  # Automatically determine ticks for the second axis
        )
    )

    # Set opacity for the traces
    fig.update_traces(opacity=0.75)

    return fig

def plot_distribution_absolute(sequences, limit):
    fig = go.Figure()
    for p in sequences:
        title = p['type']
        if 'value_counts_abs' in p:
            df = pd.DataFrame(p['value_counts_abs']).head(limit)
            n_occurrences = np.round(np.sum([item['total_occurrences'] for item in p['peaks']]),2) if len(p['peaks']) > 0 else 0
            peaks_info = ""
            for item in p['peaks']:
                idx = item['peak_index']
                occurrences = item['total_occurrences']
                peaks_info += f'{idx}({occurrences})  '
            fig.add_trace(go.Bar(
                x=df['index'],
                y=df['n_occurrences'],
                name=f'{title}: Total Amount = {n_occurrences}, PEAKS: {peaks_info}',
            )
            )
    fig.update_layout(
        width=1000,
        height=600,
        legend=dict(
            x=0,
            y=1.0,
            xanchor='left',
            yanchor='bottom'
        ),
        barmode='overlay',
        title={
            'text': "Positional distribution of occurrences in absolute numbers",
            'y': 0.05,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'bottom'}
    )
    fig.update_layout(barmode='overlay')
    fig.update_traces(opacity=0.75)
    return fig

def make_peaks_subplots(sequences):
    specs = [[{"type": "table"}] for i in range(0, len(sequences))]
    titles = []
    for item in sequences:
        seq_type = item['type']
        if 'average_peaks_distancev' in item:
            dist = item['average_peaks_distance']
            titles.append(f'{seq_type}: average distance between peaks is {round(dist)}')
        else:
            titles.append(seq_type)
    fig = make_subplots(
        rows=len(sequences), cols=1,
        shared_xaxes=True,
        vertical_spacing=0.03,
        specs=specs,
        subplot_titles=titles
        # subplot_titles=[i['type'] for i in sequences]
    )

    for i, item in enumerate(sequences, start=1):
        sequence_type = item['type']
        peaks_df = pd.DataFrame(item['peaks'])
        selected_columns = ['peak_index', 'left_bases', 'right_bases', 'total_proportion', 'total_occurrences',
                            'peak_dist']
        # Get the existing columns in peaks_df that are also in selected_columns
        existing_columns = [col for col in peaks_df.columns if col in selected_columns]

        peaks_df = peaks_df[existing_columns]

        fig.add_trace(
            go.Table(
            header=dict(values=peaks_df.columns),
            cells=dict(values=peaks_df.transpose().values)
        ),
            row=i, col=1
        )

    fig.update_layout(width=1000, height=1000, title_text='Subplots')
    return fig