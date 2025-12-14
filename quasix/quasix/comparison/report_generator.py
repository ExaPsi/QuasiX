"""
HTML report generation for CD vs AC comparison results.

This module provides tools for generating publication-quality HTML reports
with embedded plots and statistics.
"""

import os
import base64
from io import BytesIO
from typing import Dict, List, Optional, Any
from datetime import datetime
import json

import numpy as np

# Try to import matplotlib for plot embedding
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from .cd_vs_ac import ComparisonResult, ComparisonConfig


class HTMLReportGenerator:
    """Generate HTML reports for CD vs AC comparisons"""

    def __init__(self, config: ComparisonConfig = None):
        """Initialize report generator"""
        self.config = config or ComparisonConfig()
        self.sections = []
        self.plots = {}
        self.summary = {}

    def add_summary(self, summary: Dict):
        """Add summary section to report"""
        self.summary = summary

    def add_molecule_section(self, result: ComparisonResult):
        """Add a molecule result section"""
        self.sections.append(result)

    def add_plots(self, plots: Dict):
        """Add plots to report"""
        self.plots = plots

    def _figure_to_base64(self, fig) -> str:
        """Convert matplotlib figure to base64 string"""
        if not HAS_MATPLOTLIB:
            return ""

        buf = BytesIO()
        fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode('utf-8')
        plt.close(fig)
        return f"data:image/png;base64,{img_base64}"

    def _generate_css(self) -> str:
        """Generate CSS styles for the report"""
        return """
        <style>
            * {
                margin: 0;
                padding: 0;
                box-sizing: border-box;
            }

            body {
                font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto,
                           'Helvetica Neue', Arial, sans-serif;
                line-height: 1.6;
                color: #333;
                background: #f5f5f5;
                padding: 20px;
            }

            .container {
                max-width: 1200px;
                margin: 0 auto;
                background: white;
                border-radius: 10px;
                box-shadow: 0 2px 10px rgba(0,0,0,0.1);
                padding: 30px;
            }

            h1 {
                color: #2c3e50;
                border-bottom: 3px solid #3498db;
                padding-bottom: 10px;
                margin-bottom: 30px;
            }

            h2 {
                color: #34495e;
                margin-top: 30px;
                margin-bottom: 20px;
                padding-left: 10px;
                border-left: 4px solid #3498db;
            }

            h3 {
                color: #34495e;
                margin-top: 20px;
                margin-bottom: 15px;
            }

            .summary-box {
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 25px;
                border-radius: 10px;
                margin-bottom: 30px;
            }

            .summary-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 20px;
                margin-top: 20px;
            }

            .summary-item {
                background: rgba(255, 255, 255, 0.1);
                padding: 15px;
                border-radius: 8px;
                backdrop-filter: blur(10px);
            }

            .summary-label {
                font-size: 0.9em;
                opacity: 0.9;
                margin-bottom: 5px;
            }

            .summary-value {
                font-size: 1.8em;
                font-weight: bold;
            }

            .status-pass {
                color: #2ecc71;
            }

            .status-warning {
                color: #f39c12;
            }

            .status-fail {
                color: #e74c3c;
            }

            table {
                width: 100%;
                border-collapse: collapse;
                margin: 20px 0;
                background: white;
                box-shadow: 0 1px 3px rgba(0,0,0,0.1);
            }

            thead {
                background: #34495e;
                color: white;
            }

            th, td {
                padding: 12px 15px;
                text-align: left;
                border-bottom: 1px solid #ddd;
            }

            tbody tr:hover {
                background: #f8f9fa;
            }

            tbody tr:nth-child(even) {
                background: #f8f9fa;
            }

            .molecule-section {
                margin: 30px 0;
                padding: 20px;
                background: #fcfcfc;
                border: 1px solid #e0e0e0;
                border-radius: 8px;
            }

            .metrics-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
                gap: 15px;
                margin: 20px 0;
            }

            .metric-box {
                padding: 12px;
                background: white;
                border: 1px solid #e0e0e0;
                border-radius: 6px;
                text-align: center;
            }

            .metric-label {
                font-size: 0.85em;
                color: #666;
                margin-bottom: 5px;
            }

            .metric-value {
                font-size: 1.3em;
                font-weight: 600;
                color: #2c3e50;
            }

            .plot-container {
                margin: 20px 0;
                text-align: center;
            }

            .plot-container img {
                max-width: 100%;
                height: auto;
                border-radius: 8px;
                box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            }

            .plot-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
                gap: 20px;
                margin: 20px 0;
            }

            .alert {
                padding: 15px;
                margin: 20px 0;
                border-radius: 6px;
                border-left: 4px solid;
            }

            .alert-success {
                background: #d4edda;
                border-color: #28a745;
                color: #155724;
            }

            .alert-warning {
                background: #fff3cd;
                border-color: #ffc107;
                color: #856404;
            }

            .alert-danger {
                background: #f8d7da;
                border-color: #dc3545;
                color: #721c24;
            }

            .footer {
                margin-top: 50px;
                padding-top: 20px;
                border-top: 1px solid #e0e0e0;
                text-align: center;
                color: #666;
                font-size: 0.9em;
            }

            .outlier-list {
                list-style: none;
                padding: 10px;
                background: #fff5f5;
                border: 1px solid #ffcccc;
                border-radius: 6px;
                margin: 10px 0;
            }

            .outlier-item {
                padding: 5px 0;
                color: #cc0000;
            }

            .convergence-info {
                display: flex;
                justify-content: space-around;
                margin: 15px 0;
            }

            .convergence-box {
                flex: 1;
                margin: 0 10px;
                padding: 10px;
                background: #f0f0f0;
                border-radius: 6px;
                text-align: center;
            }

            code {
                background: #f4f4f4;
                padding: 2px 6px;
                border-radius: 3px;
                font-family: 'Courier New', monospace;
            }

            @media print {
                body {
                    background: white;
                }
                .container {
                    box-shadow: none;
                }
            }
        </style>
        """

    def _generate_javascript(self) -> str:
        """Generate JavaScript for interactive features"""
        return """
        <script>
            // Add sorting to tables
            function sortTable(tableId, columnIndex) {
                const table = document.getElementById(tableId);
                const tbody = table.querySelector('tbody');
                const rows = Array.from(tbody.querySelectorAll('tr'));

                const isNumeric = !isNaN(rows[0].cells[columnIndex].textContent);

                rows.sort((a, b) => {
                    const aVal = a.cells[columnIndex].textContent;
                    const bVal = b.cells[columnIndex].textContent;

                    if (isNumeric) {
                        return parseFloat(aVal) - parseFloat(bVal);
                    }
                    return aVal.localeCompare(bVal);
                });

                tbody.innerHTML = '';
                rows.forEach(row => tbody.appendChild(row));
            }

            // Highlight rows on hover
            document.addEventListener('DOMContentLoaded', function() {
                const tables = document.querySelectorAll('table');
                tables.forEach(table => {
                    table.querySelectorAll('th').forEach((th, index) => {
                        th.style.cursor = 'pointer';
                        th.onclick = () => sortTable(table.id, index);
                    });
                });
            });

            // Copy data to clipboard
            function copyToClipboard(text) {
                navigator.clipboard.writeText(text).then(() => {
                    alert('Data copied to clipboard!');
                });
            }
        </script>
        """

    def _format_summary_section(self) -> str:
        """Format the summary section"""
        if not self.summary:
            return ""

        # Determine overall status
        pass_rate = self.summary.get('pass_rate', 0)
        if pass_rate >= 0.95:
            status_class = "status-pass"
            status_text = "PASSED"
        elif pass_rate >= 0.8:
            status_class = "status-warning"
            status_text = "PARTIAL"
        else:
            status_class = "status-fail"
            status_text = "FAILED"

        html = f"""
        <div class="summary-box">
            <h2 style="color: white; border: none;">Summary Statistics</h2>
            <div class="summary-grid">
                <div class="summary-item">
                    <div class="summary-label">Overall Status</div>
                    <div class="summary-value {status_class}">{status_text}</div>
                </div>
                <div class="summary-item">
                    <div class="summary-label">Overall MAD</div>
                    <div class="summary-value">{self.summary.get('overall_mad', 0):.3f} eV</div>
                </div>
                <div class="summary-item">
                    <div class="summary-label">Pass Rate</div>
                    <div class="summary-value">{pass_rate*100:.1f}%</div>
                </div>
                <div class="summary-item">
                    <div class="summary-label">Avg Speedup</div>
                    <div class="summary-value">{self.summary.get('avg_speedup', 1):.2f}×</div>
                </div>
            </div>
        </div>
        """

        # Add detailed statistics table
        html += """
        <h3>Statistical Summary</h3>
        <table id="summary-table">
            <thead>
                <tr>
                    <th>Metric</th>
                    <th>Value</th>
                    <th>Unit</th>
                </tr>
            </thead>
            <tbody>
        """

        metrics = [
            ('Number of Molecules', self.summary.get('n_molecules', 0), ''),
            ('Overall MAD', f"{self.summary.get('overall_mad', 0):.4f}", 'eV'),
            ('Overall RMSD', f"{self.summary.get('overall_rmsd', 0):.4f}", 'eV'),
            ('MAD Threshold', f"{self.summary.get('mad_threshold', 0.05):.3f}", 'eV'),
            ('Max MAD', f"{self.summary.get('max_mad', 0):.4f}", 'eV'),
            ('Min MAD', f"{self.summary.get('min_mad', 0):.4f}", 'eV'),
            ('MAD Std Dev', f"{self.summary.get('mad_std', 0):.4f}", 'eV'),
        ]

        for label, value, unit in metrics:
            html += f"""
                <tr>
                    <td>{label}</td>
                    <td>{value}</td>
                    <td>{unit}</td>
                </tr>
            """

        html += """
            </tbody>
        </table>
        """

        return html

    def _format_molecule_section(self, result: ComparisonResult) -> str:
        """Format a single molecule result section"""
        # Determine status
        mad_threshold = self.config.mad_threshold
        if result.mad < mad_threshold:
            status_icon = "✓"
            alert_class = "alert-success"
            status_text = "Methods agree within threshold"
        elif result.mad < mad_threshold * 2:
            status_icon = "⚠"
            alert_class = "alert-warning"
            status_text = "Methods show moderate deviation"
        else:
            status_icon = "✗"
            alert_class = "alert-danger"
            status_text = "Significant deviation detected"

        html = f"""
        <div class="molecule-section">
            <h3>{status_icon} {result.molecule_name} / {result.basis_set}</h3>

            <div class="alert {alert_class}">
                {status_text}: MAD = {result.mad:.3f} eV
            </div>

            <div class="metrics-grid">
                <div class="metric-box">
                    <div class="metric-label">MAD</div>
                    <div class="metric-value">{result.mad:.3f} eV</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">RMSD</div>
                    <div class="metric-value">{result.rmsd:.3f} eV</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">Max Dev</div>
                    <div class="metric-value">{result.max_deviation:.3f} eV</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">R²</div>
                    <div class="metric-value">{result.r_squared:.4f}</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">CD Time</div>
                    <div class="metric-value">{result.cd_timing:.2f} s</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">AC Time</div>
                    <div class="metric-value">{result.ac_timing:.2f} s</div>
                </div>
            </div>

            <div class="convergence-info">
                <div class="convergence-box">
                    <strong>CD Convergence</strong><br>
                    Status: {"✓ Converged" if result.cd_converged else "✗ Not converged"}<br>
                    Iterations: {result.cd_iterations}<br>
                    Memory: {result.cd_memory_peak:.1f} MB
                </div>
                <div class="convergence-box">
                    <strong>AC Convergence</strong><br>
                    Status: {"✓ Converged" if result.ac_converged else "✗ Not converged"}<br>
                    Iterations: {result.ac_iterations}<br>
                    Memory: {result.ac_memory_peak:.1f} MB<br>
                    Condition: {result.ac_condition_number:.1f}
                </div>
            </div>
        """

        # Add confidence interval
        html += f"""
            <p><strong>95% Confidence Interval:</strong>
            MAD = {result.mad:.3f} [{result.mad_ci_lower:.3f}, {result.mad_ci_upper:.3f}] eV</p>
        """

        # Add outliers if present
        if result.outlier_indices:
            html += """
                <h4>Outliers Detected</h4>
                <ul class="outlier-list">
            """
            for idx, orbital in zip(result.outlier_indices, result.outlier_orbitals):
                cd_e = result.cd_qp_energies[idx] * 27.211  # Convert to eV
                ac_e = result.ac_qp_energies[idx] * 27.211
                diff = cd_e - ac_e
                html += f"""
                    <li class="outlier-item">
                        {orbital}: CD = {cd_e:.3f} eV, AC = {ac_e:.3f} eV, Δ = {diff:.3f} eV
                    </li>
                """
            html += "</ul>"

        # Add Z-factor comparison
        html += """
            <h4>Z-Factor Statistics</h4>
            <table>
                <thead>
                    <tr>
                        <th>Method</th>
                        <th>Mean Z</th>
                        <th>Min Z</th>
                        <th>Max Z</th>
                        <th>Std Dev</th>
                    </tr>
                </thead>
                <tbody>
        """

        for method, z_factors in [('CD', result.cd_z_factors), ('AC', result.ac_z_factors)]:
            html += f"""
                    <tr>
                        <td>{method}</td>
                        <td>{np.mean(z_factors):.3f}</td>
                        <td>{np.min(z_factors):.3f}</td>
                        <td>{np.max(z_factors):.3f}</td>
                        <td>{np.std(z_factors):.3f}</td>
                    </tr>
            """

        html += """
                </tbody>
            </table>
        </div>
        """

        return html

    def _format_plots_section(self) -> str:
        """Format the plots section"""
        if not self.plots:
            return ""

        html = "<h2>Visualization</h2>"
        html += '<div class="plot-grid">'

        for plot_name, fig in self.plots.items():
            if fig is not None:
                img_data = self._figure_to_base64(fig)
                if img_data:
                    html += f"""
                    <div class="plot-container">
                        <h3>{plot_name.replace('_', ' ').title()}</h3>
                        <img src="{img_data}" alt="{plot_name}">
                    </div>
                    """

        html += "</div>"
        return html

    def _format_results_table(self) -> str:
        """Format the results comparison table"""
        if not self.sections:
            return ""

        html = """
        <h2>Detailed Results</h2>
        <table id="results-table">
            <thead>
                <tr>
                    <th>Molecule</th>
                    <th>Basis</th>
                    <th>MAD (eV)</th>
                    <th>RMSD (eV)</th>
                    <th>R²</th>
                    <th>CD Time (s)</th>
                    <th>AC Time (s)</th>
                    <th>Speedup</th>
                    <th>Status</th>
                </tr>
            </thead>
            <tbody>
        """

        for result in self.sections:
            speedup = result.ac_timing / result.cd_timing if result.cd_timing > 0 else 0
            status = "✓" if result.mad < self.config.mad_threshold else "✗"
            status_class = "status-pass" if result.mad < self.config.mad_threshold else "status-fail"

            html += f"""
                <tr>
                    <td>{result.molecule_name}</td>
                    <td>{result.basis_set}</td>
                    <td>{result.mad:.3f}</td>
                    <td>{result.rmsd:.3f}</td>
                    <td>{result.r_squared:.4f}</td>
                    <td>{result.cd_timing:.2f}</td>
                    <td>{result.ac_timing:.2f}</td>
                    <td>{speedup:.2f}×</td>
                    <td class="{status_class}">{status}</td>
                </tr>
            """

        html += """
            </tbody>
        </table>
        """

        return html

    def generate_html(self) -> str:
        """Generate complete HTML report"""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        html = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>CD vs AC Comparison Report - {timestamp}</title>
            {self._generate_css()}
            {self._generate_javascript()}
        </head>
        <body>
            <div class="container">
                <h1>CD vs AC Method Comparison Report</h1>
                <p style="color: #666; margin-bottom: 20px;">
                    Generated on {timestamp} | QuasiX GW/BSE Implementation
                </p>

                {self._format_summary_section()}
                {self._format_results_table()}
                {self._format_plots_section()}

                <h2>Individual Molecule Analysis</h2>
                {"".join(self._format_molecule_section(r) for r in self.sections)}

                <div class="footer">
                    <p>Report generated by QuasiX Comparison Module</p>
                    <p>Contour Deformation vs Analytic Continuation GW Methods</p>
                </div>
            </div>
        </body>
        </html>
        """

        return html

    def save(self, output_path: str):
        """Save HTML report to file"""
        os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)

        html_content = self.generate_html()

        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)

        print(f"HTML report saved to: {output_path}")