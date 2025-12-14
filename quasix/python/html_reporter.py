#!/usr/bin/env python3
"""HTML report generator for GW100 benchmark results.

This module generates interactive HTML reports with charts and tables
for benchmark validation results.
"""

import json
import base64
import io
from pathlib import Path
from typing import List, Dict, Optional, Any
from datetime import datetime
import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_PLOTTING = True
except ImportError:
    HAS_PLOTTING = False

from .benchmarks import MoleculeResult, BenchmarkStatistics


class HTMLReporter:
    """Generate interactive HTML reports for benchmark results."""

    def __init__(self, title: str = "GW100 Mini-Benchmark Report"):
        self.title = title
        self.timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def generate(self,
                 results: List[MoleculeResult],
                 statistics: BenchmarkStatistics,
                 validation: Dict,
                 config: Dict,
                 output_file: Optional[Path] = None) -> str:
        """Generate complete HTML report.

        Args:
            results: List of molecule benchmark results
            statistics: Aggregated statistics
            validation: Validation results dictionary
            config: Benchmark configuration
            output_file: Optional path to save HTML file

        Returns:
            HTML string
        """
        html = self._generate_html(results, statistics, validation, config)

        if output_file:
            output_file = Path(output_file)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            with open(output_file, 'w') as f:
                f.write(html)
            print(f"HTML report saved to: {output_file}")

        return html

    def _generate_html(self, results: List[MoleculeResult],
                      statistics: BenchmarkStatistics,
                      validation: Dict,
                      config: Dict) -> str:
        """Generate the HTML content."""

        # Generate plots as base64
        plots = self._generate_plots(results, statistics) if HAS_PLOTTING else {}

        # Build HTML
        html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{self.title}</title>
    {self._get_styles()}
    {self._get_scripts()}
</head>
<body>
    <div class="container">
        <header>
            <h1>{self.title}</h1>
            <p class="timestamp">Generated: {self.timestamp}</p>
        </header>

        {self._generate_summary_section(statistics, validation)}
        {self._generate_validation_section(validation, config)}
        {self._generate_results_table(results, config)}
        {self._generate_plots_section(plots)}
        {self._generate_convergence_section(results)}
        {self._generate_metadata_section(config)}
    </div>
</body>
</html>
"""
        return html

    def _get_styles(self) -> str:
        """Get CSS styles."""
        return """
    <style>
        :root {
            --primary-color: #2c3e50;
            --secondary-color: #3498db;
            --success-color: #27ae60;
            --danger-color: #e74c3c;
            --warning-color: #f39c12;
            --bg-color: #ecf0f1;
            --card-bg: white;
            --text-color: #2c3e50;
            --border-color: #bdc3c7;
        }

        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background-color: var(--bg-color);
            color: var(--text-color);
            margin: 0;
            padding: 0;
            line-height: 1.6;
        }

        .container {
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
        }

        header {
            background: linear-gradient(135deg, var(--primary-color), var(--secondary-color));
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }

        h1 {
            margin: 0;
            font-size: 2.5em;
        }

        .timestamp {
            opacity: 0.9;
            margin-top: 10px;
        }

        .card {
            background: var(--card-bg);
            padding: 25px;
            border-radius: 8px;
            margin-bottom: 25px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }

        .card h2 {
            margin-top: 0;
            color: var(--primary-color);
            border-bottom: 2px solid var(--secondary-color);
            padding-bottom: 10px;
        }

        .summary-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }

        .metric-card {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
        }

        .metric-value {
            font-size: 2em;
            font-weight: bold;
            margin: 10px 0;
        }

        .metric-label {
            font-size: 0.9em;
            opacity: 0.9;
        }

        .status-badge {
            display: inline-block;
            padding: 5px 15px;
            border-radius: 20px;
            font-weight: bold;
            margin: 5px;
        }

        .status-pass {
            background-color: var(--success-color);
            color: white;
        }

        .status-fail {
            background-color: var(--danger-color);
            color: white;
        }

        .status-warning {
            background-color: var(--warning-color);
            color: white;
        }

        table {
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
        }

        th, td {
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid var(--border-color);
        }

        th {
            background-color: var(--primary-color);
            color: white;
            font-weight: 600;
        }

        tr:hover {
            background-color: #f5f5f5;
        }

        .plot-container {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(500px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }

        .plot-image {
            width: 100%;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }

        .validation-item {
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 10px;
            margin: 10px 0;
            background: #f8f9fa;
            border-radius: 5px;
        }

        .validation-item.pass {
            border-left: 4px solid var(--success-color);
        }

        .validation-item.fail {
            border-left: 4px solid var(--danger-color);
        }

        .progress-bar {
            width: 100%;
            height: 30px;
            background-color: #e0e0e0;
            border-radius: 15px;
            overflow: hidden;
            margin: 20px 0;
        }

        .progress-fill {
            height: 100%;
            background: linear-gradient(90deg, var(--success-color), var(--secondary-color));
            display: flex;
            align-items: center;
            justify-content: center;
            color: white;
            font-weight: bold;
            transition: width 0.3s ease;
        }
    </style>
        """

    def _get_scripts(self) -> str:
        """Get JavaScript for interactivity."""
        return """
    <script>
        // Add sorting functionality to tables
        document.addEventListener('DOMContentLoaded', function() {
            const tables = document.querySelectorAll('table');
            tables.forEach(table => {
                const headers = table.querySelectorAll('th');
                headers.forEach((header, index) => {
                    header.style.cursor = 'pointer';
                    header.addEventListener('click', () => sortTable(table, index));
                });
            });
        });

        function sortTable(table, column) {
            const tbody = table.querySelector('tbody');
            const rows = Array.from(tbody.querySelectorAll('tr'));
            const isNumeric = !isNaN(rows[0].cells[column].textContent);

            rows.sort((a, b) => {
                const aVal = a.cells[column].textContent;
                const bVal = b.cells[column].textContent;

                if (isNumeric) {
                    return parseFloat(aVal) - parseFloat(bVal);
                } else {
                    return aVal.localeCompare(bVal);
                }
            });

            // Toggle sort direction
            if (table.dataset.sortColumn == column) {
                rows.reverse();
                table.dataset.sortColumn = -1;
            } else {
                table.dataset.sortColumn = column;
            }

            // Re-append rows
            rows.forEach(row => tbody.appendChild(row));
        }
    </script>
        """

    def _generate_summary_section(self, statistics: BenchmarkStatistics,
                                  validation: Dict) -> str:
        """Generate summary section."""
        status_class = "status-pass" if validation['overall_pass'] else "status-fail"
        status_text = "PASSED" if validation['overall_pass'] else "FAILED"

        return f"""
        <div class="card">
            <h2>Summary</h2>
            <div style="text-align: center; margin: 20px 0;">
                <span class="status-badge {status_class}" style="font-size: 1.5em;">
                    Overall Status: {status_text}
                </span>
            </div>

            <div class="summary-grid">
                <div class="metric-card">
                    <div class="metric-label">Molecules Tested</div>
                    <div class="metric-value">{statistics.n_molecules}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Mean Absolute Deviation</div>
                    <div class="metric-value">{statistics.mad:.3f} eV</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">RMSD</div>
                    <div class="metric-value">{statistics.rmsd:.3f} eV</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">R²</div>
                    <div class="metric-value">{statistics.r_squared:.4f}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Pass Rate</div>
                    <div class="metric-value">{statistics.pass_rate*100:.1f}%</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Max Deviation</div>
                    <div class="metric-value">{statistics.max_dev:.3f} eV</div>
                </div>
            </div>

            <div class="progress-bar">
                <div class="progress-fill" style="width: {statistics.pass_rate*100}%;">
                    {statistics.pass_rate*100:.1f}% Tests Passed
                </div>
            </div>
        </div>
        """

    def _generate_validation_section(self, validation: Dict, config: Dict) -> str:
        """Generate validation criteria section."""
        items = []

        criteria = [
            ("MAD ≤ {:.2f} eV".format(config.get('max_deviation_ev', 0.2)),
             validation['mad_pass']),
            ("RMSD ≤ {:.2f} eV".format(config.get('max_deviation_ev', 0.2) * 1.2),
             validation['rmsd_pass']),
            ("R² ≥ 0.95", validation['correlation_pass']),
            ("Pass Rate ≥ 75%", validation['pass_rate_pass'])
        ]

        for criterion, passed in criteria:
            status = "pass" if passed else "fail"
            icon = "✓" if passed else "✗"
            items.append(f"""
                <div class="validation-item {status}">
                    <span>{criterion}</span>
                    <span style="font-size: 1.5em;">{icon}</span>
                </div>
            """)

        return f"""
        <div class="card">
            <h2>Validation Criteria</h2>
            {''.join(items)}
        </div>
        """

    def _generate_results_table(self, results: List[MoleculeResult],
                                config: Dict) -> str:
        """Generate results table."""
        rows = []
        for r in results:
            passed = r.passes_threshold(config.get('max_deviation_ev', 0.2))
            status_icon = "✓" if passed else "✗"
            row_class = "" if passed else 'style="background-color: #ffe6e6;"'

            ea_str = f"{r.ea_calc:.3f}" if r.ea_calc is not None else "N/A"
            ea_ref_str = f"{r.ea_ref:.3f}" if r.ea_ref is not None else "N/A"

            rows.append(f"""
                <tr {row_class}>
                    <td>{r.molecule}</td>
                    <td>{r.ip_calc:.3f}</td>
                    <td>{r.ip_ref:.3f}</td>
                    <td>{r.ip_deviation:.3f}</td>
                    <td>{ea_str}</td>
                    <td>{ea_ref_str}</td>
                    <td>{r.elapsed_seconds:.2f}</td>
                    <td style="text-align: center;">{status_icon}</td>
                </tr>
            """)

        return f"""
        <div class="card">
            <h2>Individual Results</h2>
            <table>
                <thead>
                    <tr>
                        <th>Molecule</th>
                        <th>IP Calc (eV)</th>
                        <th>IP Ref (eV)</th>
                        <th>IP Dev (eV)</th>
                        <th>EA Calc (eV)</th>
                        <th>EA Ref (eV)</th>
                        <th>Time (s)</th>
                        <th>Pass</th>
                    </tr>
                </thead>
                <tbody>
                    {''.join(rows)}
                </tbody>
            </table>
        </div>
        """

    def _generate_plots(self, results: List[MoleculeResult],
                       statistics: BenchmarkStatistics) -> Dict[str, str]:
        """Generate plots as base64 encoded images."""
        if not HAS_PLOTTING:
            return {}

        plots = {}

        # Set style
        sns.set_style("whitegrid")

        # 1. Correlation plot
        fig, ax = plt.subplots(figsize=(8, 6))
        ip_ref = np.array([r.ip_ref for r in results])
        ip_calc = np.array([r.ip_calc for r in results])

        ax.scatter(ip_ref, ip_calc, s=100, alpha=0.7, edgecolors='black')

        # Fit line
        z = np.polyfit(ip_ref, ip_calc, 1)
        p = np.poly1d(z)
        x_line = np.linspace(ip_ref.min(), ip_ref.max(), 100)
        ax.plot(x_line, p(x_line), 'r--', alpha=0.8,
               label=f'y = {z[0]:.3f}x + {z[1]:.3f}')
        ax.plot(x_line, x_line, 'k:', alpha=0.5, label='y = x')

        # Annotate points
        for r in results:
            ax.annotate(r.molecule, (r.ip_ref, r.ip_calc),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, alpha=0.7)

        ax.set_xlabel('Reference IP (eV)', fontsize=12)
        ax.set_ylabel('Calculated IP (eV)', fontsize=12)
        ax.set_title('Ionization Potential Correlation', fontsize=14)
        ax.legend()
        ax.grid(True, alpha=0.3)

        plots['correlation'] = self._fig_to_base64(fig)
        plt.close()

        # 2. Deviation bar plot
        fig, ax = plt.subplots(figsize=(10, 6))
        molecules = [r.molecule for r in results]
        deviations = [r.ip_deviation for r in results]
        colors = ['green' if abs(d) <= 0.2 else 'red' for d in deviations]

        bars = ax.bar(range(len(molecules)), deviations, color=colors, alpha=0.7)
        ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
        ax.axhline(y=0.2, color='r', linestyle='--', alpha=0.5)
        ax.axhline(y=-0.2, color='r', linestyle='--', alpha=0.5)

        ax.set_xticks(range(len(molecules)))
        ax.set_xticklabels(molecules, rotation=45, ha='right')
        ax.set_ylabel('IP Deviation (eV)', fontsize=12)
        ax.set_title('Individual Deviations from Reference', fontsize=14)
        ax.grid(True, alpha=0.3, axis='y')

        plots['deviations'] = self._fig_to_base64(fig)
        plt.close()

        return plots

    def _fig_to_base64(self, fig) -> str:
        """Convert matplotlib figure to base64 string."""
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode('utf-8')
        buf.close()
        return f"data:image/png;base64,{img_base64}"

    def _generate_plots_section(self, plots: Dict[str, str]) -> str:
        """Generate plots section."""
        if not plots:
            return ""

        plot_html = []
        for name, data in plots.items():
            plot_html.append(f'<img src="{data}" class="plot-image" alt="{name}">')

        return f"""
        <div class="card">
            <h2>Visualizations</h2>
            <div class="plot-container">
                {''.join(plot_html)}
            </div>
        </div>
        """

    def _generate_convergence_section(self, results: List[MoleculeResult]) -> str:
        """Generate convergence analysis section."""
        conv_data = []
        for r in results:
            if r.convergence_iterations > 0:
                conv_data.append(f"""
                    <tr>
                        <td>{r.molecule}</td>
                        <td>{r.convergence_iterations}</td>
                        <td>{r.max_error:.2e}</td>
                        <td>{np.mean(r.z_factors):.3f}</td>
                        <td>{r.elapsed_seconds:.2f}</td>
                    </tr>
                """)

        if not conv_data:
            return ""

        return f"""
        <div class="card">
            <h2>Convergence Analysis</h2>
            <table>
                <thead>
                    <tr>
                        <th>Molecule</th>
                        <th>Iterations</th>
                        <th>Final Error</th>
                        <th>Mean Z-factor</th>
                        <th>Time (s)</th>
                    </tr>
                </thead>
                <tbody>
                    {''.join(conv_data)}
                </tbody>
            </table>
        </div>
        """

    def _generate_metadata_section(self, config: Dict) -> str:
        """Generate metadata section."""
        items = []
        for key, value in config.items():
            items.append(f"<tr><td><strong>{key}</strong></td><td>{value}</td></tr>")

        return f"""
        <div class="card">
            <h2>Configuration</h2>
            <table style="width: auto;">
                <tbody>
                    {''.join(items)}
                </tbody>
            </table>
        </div>
        """


def generate_html_report(results: List[MoleculeResult],
                        statistics: BenchmarkStatistics,
                        validation: Dict,
                        config: Dict,
                        output_file: Optional[Path] = None) -> str:
    """Convenience function to generate HTML report.

    Args:
        results: List of molecule benchmark results
        statistics: Aggregated statistics
        validation: Validation results
        config: Benchmark configuration
        output_file: Optional output file path

    Returns:
        HTML string
    """
    reporter = HTMLReporter()
    return reporter.generate(results, statistics, validation, config, output_file)