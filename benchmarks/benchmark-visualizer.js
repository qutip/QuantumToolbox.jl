/**
 * Benchmark Visualizer for QuantumToolbox.jl
 * 
 * This script visualizes benchmark data using Chart.js. It organizes benchmarks
 * by topic and case, allowing for easy comparison across commits.
 */
"use strict";

// Configuration for chart colors
const CONFIG = {
  colors: ["#a270ba", "#00add8", "#f1e05a", "#3572a5", "#f34b7d", "#b07219", 
           "#178600", "#38ff38", "#ff3838", "#dea584", "#000080", "#333333"],
};

class BenchmarkVisualizer {
  constructor(data) {
    this.data = data;
    this.dataSet = this.collectBenchmarks(data.entries["Benchmark Results"]);
  }

  /**
   * Initialize the visualizer: render header, charts, and footer
   */
  init() {
    this.renderHeader();
    this.renderCharts();
    this.setupFooter();
  }

  /**
   * Render the page header with last update time and repository link
   */
  renderHeader() {
    document.getElementById("last-update").textContent = new Date(this.data.lastUpdate).toString();
    const repoLink = document.getElementById("repository-link");
    repoLink.href = this.data.repoUrl;
    repoLink.textContent = this.data.repoUrl;
  }

  /**
   * Setup the footer with download button for raw JSON data
   */
  setupFooter() {
    document.getElementById("dl-button").onclick = () => {
      const link = document.createElement("a");
      link.href = "data:," + JSON.stringify(this.data, null, 2);
      link.download = "benchmark_data.json";
      link.click();
    };
  }

  /**
   * Collect and organize benchmarks into a hierarchical structure:
   * Topic -> Case -> Legend -> [Results]
   * 
   * Benchmark names are split by "/" to create this hierarchy.
   * Example: "Time Evolution/time-independent/sesolve" becomes:
   *   Topic: "Time Evolution"
   *   Case: "time-independent"
   *   Legend: "sesolve"
   * 
   * @param {Array} entries - Array of benchmark entries from data.js
   * @returns {Map} Hierarchical map of organized benchmarks
   */
  collectBenchmarks(entries) {
    const topicMap = new Map();
    
    for (const entry of entries) {
      for (const bench of entry.benches) {
        const result = { commit: entry.commit, date: entry.date, tool: entry.tool, bench };
        const parts = bench.name.split("/");
        
        // Get or create topic (e.g., "Time Evolution")
        const topic = this.getOrCreate(topicMap, parts[0], () => new Map());
        
        // Get or create case (e.g., "time-independent")
        const caseName = parts[1] || "default";
        const caseData = this.getOrCreate(topic, caseName, () => new Map());
        
        // Get or create legend (e.g., "sesolve")
        const legendName = parts.slice(2).join("/") || "default";
        const legendData = this.getOrCreate(caseData, legendName, () => []);
        
        // Add the result to the legend's data array
        legendData.push(result);
      }
    }
    
    return topicMap;
  }

  /**
   * Helper function to get or create a map entry
   * 
   * @param {Map} map - The map to search/modify
   * @param {string} key - The key to look for
   * @param {Function} creator - Function that creates the default value if key doesn't exist
   * @returns {*} The existing or newly created value
   */
  getOrCreate(map, key, creator) {
    if (!map.has(key)) map.set(key, creator());
    return map.get(key);
  }

  /**
   * Render all charts organized by topics
   */
  renderCharts() {
    const container = document.createElement("div");
    container.className = "charts-container";
    document.getElementById("main").appendChild(container);
    
    for (const [topicName, topicData] of this.dataSet.entries()) {
      this.renderTopic(topicName, topicData, container);
    }
  }

  /**
   * Render a topic section with its title and grid of charts
   * 
   * @param {string} topicName - Name of the topic (e.g., "Time Evolution")
   * @param {Map} topicData - Map of cases within this topic
   * @param {HTMLElement} container - Parent container to append to
   */
  renderTopic(topicName, topicData, container) {
    const section = document.createElement("div");
    section.className = "benchmark-set";
    
    const title = document.createElement("h1");
    title.className = "benchmark-title";
    title.textContent = topicName;
    section.appendChild(title);
    
    const grid = document.createElement("div");
    grid.className = "charts-grid";
    section.appendChild(grid);
    
    container.appendChild(section);
    
    // Render each case as a separate chart
    for (const [caseName, caseData] of topicData.entries()) {
      this.renderChart(caseName, caseData, grid);
    }
  }

  /**
   * Render a single chart for a case
   * 
   * @param {string} caseName - Name of the case (chart title)
   * @param {Map} dataset - Map of legend names to data arrays
   * @param {HTMLElement} parent - Parent element to append chart to
   */
  renderChart(caseName, dataset, parent) {
    // Create chart container
    const container = document.createElement("div");
    container.className = "chartCard";
    const canvas = document.createElement("canvas");
    canvas.className = "benchmark-chart";
    container.appendChild(canvas);
    parent.appendChild(container);
    
    // Get all unique commit IDs across all legends
    const commitIds = [...new Set([...dataset.values()].flat().map(d => d.commit.id))];
    
    // Build datasets for Chart.js (one dataset per legend/line)
    const datasets = [];
    let colorIdx = 0;
    
    for (const [legendName, legendData] of dataset.entries()) {
      const color = CONFIG.colors[colorIdx++ % CONFIG.colors.length];
      
      // Create a map for fast lookup: commitId -> benchmark value
      const dataMap = new Map(legendData.map(d => [d.commit.id, d.bench.value / 1e6]));
      
      datasets.push({
        label: legendName,
        data: commitIds.map(id => dataMap.get(id) ?? null),
        borderColor: color,
        backgroundColor: color + "40", // Add transparency to background
        fill: false,
        tension: 0.1,
        pointRadius: 3,
        pointHoverRadius: 5,
      });
    }
    
    // Create the Chart.js chart
    new Chart(canvas, {
      type: "line",
      data: { 
        labels: commitIds.map(id => id.slice(0, 7)), // Short commit hashes
        datasets 
      },
      options: {
        responsive: true,
        maintainAspectRatio: false, // Allow chart to fill container height
        plugins: {
          title: { display: true, text: caseName, font: { size: 16, weight: 'bold' } },
          tooltip: { callbacks: this.createTooltipCallbacks(dataset, commitIds) },
          legend: { display: true, position: 'top' },
        },
        scales: {
          x: { title: { text: "Commit", display: true } },
          y: { title: { text: "Time (ms)", display: true }, beginAtZero: true },
        },
        onClick: (e, items) => this.handleClick(items, commitIds, dataset),
      },
    });
  }

  /**
   * Create custom tooltip callbacks to show detailed benchmark info
   * 
   * @param {Map} dataset - The dataset for this chart
   * @param {Array} commitIds - Array of commit IDs
   * @returns {Object} Chart.js tooltip callback configuration
   */
  createTooltipCallbacks(dataset, commitIds) {
    return {
      // Show commit message and author after title
      afterTitle: items => {
        const d = this.getDataPoint(items[0], dataset, commitIds);
        return d ? `\n${d.commit.message}\n\n${d.commit.timestamp} by @${d.commit.committer.username}\n` : '';
      },
      
      // Show benchmark value with unit and optional range
      label: context => {
        const d = this.getDataPoint(context, dataset, commitIds);
        if (!d) return context.dataset.label;
        let label = `${context.dataset.label} : ${d.bench.value} ${d.bench.unit}`;
        if (d.bench.range) label += ` (${d.bench.range})`;
        return label;
      },
      
      // Show extra benchmark info if available
      afterLabel: item => {
        const d = this.getDataPoint(item, dataset, commitIds);
        return d?.bench.extra ? `\n${d.bench.extra}` : '';
      },
    };
  }

  /**
   * Get the data point for a specific chart element
   * 
   * @param {Object} item - Chart.js item (from tooltip or click event)
   * @param {Map} dataset - The dataset for this chart
   * @param {Array} commitIds - Array of commit IDs
   * @returns {Object|undefined} The data point or undefined if not found
   */
  getDataPoint(item, dataset, commitIds) {
    const label = item.dataset?.label || item.element?.$context?.dataset?.label;
    const commitId = commitIds[item.dataIndex];
    return dataset.get(label)?.find(d => d.commit.id === commitId);
  }

  /**
   * Handle chart click events - open commit URL in new tab
   * 
   * @param {Array} items - Clicked chart items
   * @param {Array} commitIds - Array of commit IDs
   * @param {Map} dataset - The dataset for this chart
   */
  handleClick(items, commitIds, dataset) {
    if (items.length === 0) return;
    const d = this.getDataPoint(items[0], dataset, commitIds);
    if (d?.commit.url) window.open(d.commit.url, '_blank');
  }
}

// Initialize visualizer when DOM is ready
document.addEventListener('DOMContentLoaded', () => {
  if (window.BENCHMARK_DATA) {
    new BenchmarkVisualizer(window.BENCHMARK_DATA).init();
  }
});
