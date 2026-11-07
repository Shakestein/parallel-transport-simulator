#pragma once

#include "vec3.h"
#include "surfaces.h"
#include "transport.h"
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================================
// Generates interactive HTML + Plotly.js visualization
// ============================================================================
class Visualization {
public:
    // Generate and save HTML visualization, returns filepath
    static std::string generate(const Surface& surface,
                                const Curve& curve,
                                const TransportResult& result,
                                const std::string& output_dir = "output",
                                int mesh_resolution = 60,
                                int arrow_skip = -1) {
        // Default: show ~25 arrows
        if (arrow_skip < 0) {
            arrow_skip = std::max(1, static_cast<int>(result.positions.size()) / 25);
        }

        // Generate surface mesh
        double u_min, u_max, v_min, v_max;
        surface.domain(u_min, u_max, v_min, v_max);

        int nu = mesh_resolution, nv = mesh_resolution;
        std::vector<std::vector<double>> X(nu, std::vector<double>(nv));
        std::vector<std::vector<double>> Y(nu, std::vector<double>(nv));
        std::vector<std::vector<double>> Z(nu, std::vector<double>(nv));

        for (int i = 0; i < nu; i++) {
            for (int j = 0; j < nv; j++) {
                double u = u_min + (u_max - u_min) * i / (nu - 1);
                double v = v_min + (v_max - v_min) * j / (nv - 1);
                Vec3 p = surface.position(u, v);
                X[i][j] = p.x;
                Y[i][j] = p.y;
                Z[i][j] = p.z;
            }
        }

        // Build HTML
        std::ostringstream html;
        html << R"(<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Parallel Transport Simulator</title>
<script src="https://cdn.plot.ly/plotly-2.35.0.min.js"></script>
<style>
  * { margin: 0; padding: 0; box-sizing: border-box; }
  body {
    font-family: 'Segoe UI', system-ui, -apple-system, sans-serif;
    background: #0a0a1a;
    color: #e0e0e0;
    overflow: hidden;
  }
  #header {
    background: linear-gradient(135deg, #1a1a3e, #0d0d2b);
    padding: 12px 24px;
    display: flex;
    justify-content: space-between;
    align-items: center;
    border-bottom: 1px solid #2a2a5a;
    box-shadow: 0 2px 12px rgba(0,0,0,0.5);
  }
  #header h1 {
    font-size: 1.3em;
    font-weight: 600;
    background: linear-gradient(90deg, #7eb8ff, #c084fc);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
  }
  #info-panel {
    display: flex;
    gap: 24px;
    font-size: 0.85em;
  }
  .info-item {
    display: flex;
    flex-direction: column;
    align-items: center;
  }
  .info-label {
    color: #888;
    font-size: 0.75em;
    text-transform: uppercase;
    letter-spacing: 1px;
  }
  .info-value {
    font-weight: 600;
    font-size: 1.1em;
    margin-top: 2px;
  }
  .info-value.good { color: #4ade80; }
  .info-value.highlight { color: #f59e0b; }
  .info-value.surface-name { color: #7eb8ff; }
  #plot { width: 100vw; height: calc(100vh - 56px); }
  #legend-overlay {
    position: absolute;
    bottom: 20px;
    left: 20px;
    background: rgba(20, 20, 50, 0.9);
    border: 1px solid #2a2a5a;
    border-radius: 8px;
    padding: 12px 16px;
    font-size: 0.8em;
    backdrop-filter: blur(10px);
  }
  .legend-item {
    display: flex;
    align-items: center;
    gap: 8px;
    margin: 4px 0;
  }
  .legend-color {
    width: 16px;
    height: 4px;
    border-radius: 2px;
  }
</style>
</head>
<body>
<div id="header">
  <h1>&#x2B50; Parallel Transport Simulator</h1>
  <div id="info-panel">
    <div class="info-item">
      <span class="info-label">Surface</span>
      <span class="info-value surface-name">)" << surface.name() << R"(</span>
    </div>
    <div class="info-item">
      <span class="info-label">Curve</span>
      <span class="info-value surface-name">)" << curve.name() << R"(</span>
    </div>
    <div class="info-item">
      <span class="info-label">Holonomy Angle</span>
      <span class="info-value highlight">)" << format_angle(result.holonomy_angle) << R"(</span>
    </div>
    <div class="info-item">
      <span class="info-label">Max Norm Drift</span>
      <span class="info-value good">)" << format_scientific(result.max_norm_drift) << R"(</span>
    </div>
    <div class="info-item">
      <span class="info-label">Steps</span>
      <span class="info-value" style="color:#c084fc">)" << result.positions.size() - 1 << R"(</span>
    </div>
  </div>
</div>
<div id="plot"></div>
<div id="legend-overlay">
  <div class="legend-item"><div class="legend-color" style="background:#ef4444"></div>Transport Path</div>
  <div class="legend-item"><div class="legend-color" style="background:#22d3ee"></div>Initial Vector</div>
  <div class="legend-item"><div class="legend-color" style="background:#fbbf24"></div>Transported Vectors</div>
  <div class="legend-item"><div class="legend-color" style="background:#f97316"></div>Final Vector</div>
</div>
<script>
)";

        // Emit surface data
        html << "var surfX = " << to_json_2d(X) << ";\n";
        html << "var surfY = " << to_json_2d(Y) << ";\n";
        html << "var surfZ = " << to_json_2d(Z) << ";\n";

        // Emit curve data
        html << "var curveX = " << to_json_1d_vec3(result.positions, 0) << ";\n";
        html << "var curveY = " << to_json_1d_vec3(result.positions, 1) << ";\n";
        html << "var curveZ = " << to_json_1d_vec3(result.positions, 2) << ";\n";

        // Emit arrow data (sampled positions and vectors)
        std::vector<Vec3> arrow_bases, arrow_dirs;
        std::vector<int> arrow_indices;
        for (size_t i = 0; i < result.positions.size(); i += arrow_skip) {
            arrow_bases.push_back(result.positions[i]);
            arrow_dirs.push_back(result.vectors_3d[i]);
            arrow_indices.push_back(static_cast<int>(i));
        }
        // Always include last point
        if (arrow_indices.empty() || arrow_indices.back() != static_cast<int>(result.positions.size() - 1)) {
            arrow_bases.push_back(result.positions.back());
            arrow_dirs.push_back(result.vectors_3d.back());
            arrow_indices.push_back(static_cast<int>(result.positions.size() - 1));
        }

        // Compute arrow scale: make arrows proportional to surface size
        double max_extent = compute_extent(result.positions);
        double arrow_scale = max_extent * 0.15;

        html << R"(
var traces = [];

// Surface mesh
traces.push({
  type: 'surface',
  x: surfX, y: surfY, z: surfZ,
  colorscale: [
    [0, 'rgba(30, 30, 80, 0.85)'],
    [0.25, 'rgba(40, 60, 120, 0.85)'],
    [0.5, 'rgba(60, 100, 160, 0.85)'],
    [0.75, 'rgba(80, 140, 180, 0.85)'],
    [1, 'rgba(120, 180, 220, 0.85)']
  ],
  opacity: 0.7,
  showscale: false,
  lighting: {
    ambient: 0.6,
    diffuse: 0.7,
    specular: 0.3,
    roughness: 0.5
  },
  contours: {
    x: { show: true, color: 'rgba(100,100,200,0.2)', width: 1 },
    y: { show: true, color: 'rgba(100,100,200,0.2)', width: 1 },
    z: { show: true, color: 'rgba(100,100,200,0.2)', width: 1 }
  },
  name: 'Surface',
  hoverinfo: 'skip'
});

// Transport curve
traces.push({
  type: 'scatter3d',
  x: curveX, y: curveY, z: curveZ,
  mode: 'lines',
  line: { color: '#ef4444', width: 7 },
  name: 'Transport Path',
  hoverinfo: 'skip'
});

)";

        // Generate arrow traces (shafts + coneheads)
        html << "// Transported vector arrows\n";
        for (size_t idx = 0; idx < arrow_bases.size(); idx++) {
            Vec3 base = arrow_bases[idx];
            Vec3 dir = arrow_dirs[idx];

            // Normalize and scale
            double dir_norm = dir.norm();
            if (dir_norm < 1e-15) continue;
            Vec3 dir_scaled = dir * (arrow_scale / dir_norm);
            Vec3 tip = base + dir_scaled;

            // Color: first=cyan, last=orange, middle=gold
            std::string color;
            if (idx == 0) {
                color = "#22d3ee";  // cyan for initial
            } else if (idx == arrow_bases.size() - 1) {
                color = "#f97316";  // orange for final
            } else {
                // Gradient from gold to amber
                double frac = static_cast<double>(idx) / (arrow_bases.size() - 1);
                int r = static_cast<int>(251 * (1.0 - frac * 0.3));
                int g = static_cast<int>(191 - frac * 60);
                int b = static_cast<int>(36 + frac * 30);
                color = rgb_to_hex(r, g, b);
            }

            double shaft_width = (idx == 0 || idx == arrow_bases.size() - 1) ? 7.0 : 5.0;

            // Arrow shaft
            html << "traces.push({\n"
                 << "  type: 'scatter3d',\n"
                 << "  x: [" << base.x << "," << tip.x << "],\n"
                 << "  y: [" << base.y << "," << tip.y << "],\n"
                 << "  z: [" << base.z << "," << tip.z << "],\n"
                 << "  mode: 'lines',\n"
                 << "  line: { color: '" << color << "', width: " << shaft_width << " },\n"
                 << "  showlegend: false,\n"
                 << "  hoverinfo: 'skip'\n"
                 << "});\n";

            // Arrow cone head
            html << "traces.push({\n"
                 << "  type: 'cone',\n"
                 << "  x: [" << tip.x << "],\n"
                 << "  y: [" << tip.y << "],\n"
                 << "  z: [" << tip.z << "],\n"
                 << "  u: [" << dir_scaled.x << "],\n"
                 << "  v: [" << dir_scaled.y << "],\n"
                 << "  w: [" << dir_scaled.z << "],\n"
                 << "  sizemode: 'absolute',\n"
                 << "  sizeref: " << arrow_scale * 0.3 << ",\n"
                 << "  colorscale: [[0, '" << color << "'], [1, '" << color << "']],\n"
                 << "  showscale: false,\n"
                 << "  showlegend: false,\n"
                 << "  hoverinfo: 'skip',\n"
                 << "  anchor: 'tip'\n"
                 << "});\n";
        }

        // Start/end marker dots
        Vec3 start_pos = result.positions.front();
        Vec3 end_pos = result.positions.back();

        html << "// Start and end markers\n";
        html << "traces.push({\n"
             << "  type: 'scatter3d',\n"
             << "  x: [" << start_pos.x << "],\n"
             << "  y: [" << start_pos.y << "],\n"
             << "  z: [" << start_pos.z << "],\n"
             << "  mode: 'markers+text',\n"
             << "  marker: { size: 7, color: '#22d3ee', symbol: 'circle' },\n"
             << "  text: ['Start'],\n"
             << "  textposition: 'top center',\n"
             << "  textfont: { color: '#22d3ee', size: 12 },\n"
             << "  showlegend: false,\n"
             << "  hoverinfo: 'text'\n"
             << "});\n";

        html << "traces.push({\n"
             << "  type: 'scatter3d',\n"
             << "  x: [" << end_pos.x << "],\n"
             << "  y: [" << end_pos.y << "],\n"
             << "  z: [" << end_pos.z << "],\n"
             << "  mode: 'markers+text',\n"
             << "  marker: { size: 7, color: '#f97316', symbol: 'diamond' },\n"
             << "  text: ['End'],\n"
             << "  textposition: 'top center',\n"
             << "  textfont: { color: '#f97316', size: 12 },\n"
             << "  showlegend: false,\n"
             << "  hoverinfo: 'text'\n"
             << "});\n";

        // Layout and render
        html << R"(
var layout = {
  paper_bgcolor: '#0a0a1a',
  plot_bgcolor: '#0a0a1a',
  margin: { l: 0, r: 0, t: 0, b: 0 },
  scene: {
    xaxis: {
      title: 'X',
      gridcolor: 'rgba(100,100,200,0.15)',
      zerolinecolor: 'rgba(100,100,200,0.3)',
      color: '#888',
      backgroundcolor: '#0a0a1a'
    },
    yaxis: {
      title: 'Y',
      gridcolor: 'rgba(100,100,200,0.15)',
      zerolinecolor: 'rgba(100,100,200,0.3)',
      color: '#888',
      backgroundcolor: '#0a0a1a'
    },
    zaxis: {
      title: 'Z',
      gridcolor: 'rgba(100,100,200,0.15)',
      zerolinecolor: 'rgba(100,100,200,0.3)',
      color: '#888',
      backgroundcolor: '#0a0a1a'
    },
    aspectmode: 'data',
    camera: {
      eye: { x: 1.5, y: 1.5, z: 1.2 },
      up: { x: 0, y: 0, z: 1 }
    }
  },
  showlegend: false
};

Plotly.newPlot('plot', traces, layout, {
  responsive: true,
  displayModeBar: true,
  modeBarButtonsToRemove: ['toImage', 'sendDataToCloud'],
  displaylogo: false
});
</script>
</body>
</html>
)";

        // Write to file
        std::string filename = sanitize_filename(surface.name()) + "_" +
                               sanitize_filename(curve.name()) + ".html";
        std::string filepath = output_dir + "/" + filename;

        // Create output directory (platform-specific)
        #ifdef _WIN32
            std::string cmd = "mkdir \"" + output_dir + "\" 2>nul";
        #else
            std::string cmd = "mkdir -p \"" + output_dir + "\"";
        #endif
        system(cmd.c_str());

        std::ofstream file(filepath);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filepath << " for writing.\n";
            return "";
        }
        file << html.str();
        file.close();

        std::cout << "  Visualization saved to: " << filepath << "\n";
        return filepath;
    }

    // Open the HTML file in the default browser
    static void open_in_browser(const std::string& filepath) {
        if (filepath.empty()) return;

        #ifdef _WIN32
            std::string cmd = "start \"\" \"" + filepath + "\"";
        #elif __APPLE__
            std::string cmd = "open \"" + filepath + "\"";
        #else
            std::string cmd = "xdg-open \"" + filepath + "\"";
        #endif
        system(cmd.c_str());
    }

private:
    // Convert 2D array to JSON
    static std::string to_json_2d(const std::vector<std::vector<double>>& arr) {
        std::ostringstream ss;
        ss.precision(8);
        ss << "[";
        for (size_t i = 0; i < arr.size(); i++) {
            ss << "[";
            for (size_t j = 0; j < arr[i].size(); j++) {
                ss << arr[i][j];
                if (j + 1 < arr[i].size()) ss << ",";
            }
            ss << "]";
            if (i + 1 < arr.size()) ss << ",";
        }
        ss << "]";
        return ss.str();
    }

    // Extract one component from Vec3 array to JSON
    static std::string to_json_1d_vec3(const std::vector<Vec3>& arr, int comp) {
        std::ostringstream ss;
        ss.precision(8);
        ss << "[";
        for (size_t i = 0; i < arr.size(); i++) {
            switch (comp) {
                case 0: ss << arr[i].x; break;
                case 1: ss << arr[i].y; break;
                case 2: ss << arr[i].z; break;
            }
            if (i + 1 < arr.size()) ss << ",";
        }
        ss << "]";
        return ss.str();
    }

    // Compute extent of point set for arrow scaling
    static double compute_extent(const std::vector<Vec3>& points) {
        if (points.empty()) return 1.0;
        double max_dist = 0.0;
        Vec3 center = {0, 0, 0};
        for (const auto& p : points) center += p;
        center = center / static_cast<double>(points.size());
        for (const auto& p : points) {
            max_dist = std::max(max_dist, (p - center).norm());
        }
        return std::max(max_dist, 0.1);
    }

    // Format angle for display
    static std::string format_angle(double rad) {
        std::ostringstream ss;
        double deg = rad * 180.0 / M_PI;
        ss.precision(4);
        ss << std::fixed << deg << "&deg; (" << std::fixed << rad << " rad)";
        return ss.str();
    }

    // Format small number in scientific notation
    static std::string format_scientific(double val) {
        std::ostringstream ss;
        ss.precision(2);
        ss << std::scientific << val;
        return ss.str();
    }

    // RGB to hex color
    static std::string rgb_to_hex(int r, int g, int b) {
        char buf[8];
        snprintf(buf, sizeof(buf), "#%02x%02x%02x", r & 0xFF, g & 0xFF, b & 0xFF);
        return std::string(buf);
    }

    // Sanitize string for filename
    static std::string sanitize_filename(const std::string& s) {
        std::string result;
        for (char c : s) {
            if (std::isalnum(c) || c == '_' || c == '-') {
                result += std::tolower(c);
            } else if (c == ' ' || c == '(' || c == ')' || c == '=' || c == ',' || c == '.') {
                if (!result.empty() && result.back() != '_')
                    result += '_';
            }
        }
        // Remove trailing underscore
        while (!result.empty() && result.back() == '_') result.pop_back();
        return result;
    }
};
