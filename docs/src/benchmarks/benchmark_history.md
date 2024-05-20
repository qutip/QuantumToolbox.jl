```@raw html
<style>
    .chart-container {
        width: 100%;
        max-width: 600px;
        margin: auto;
    }
</style>

<div class="chart-container">
    <canvas id="myChart"></canvas>
</div>

<script src="https://cdn.jsdelivr.net/npm/chart.js@2.9.2/dist/Chart.min.js"></script>
<script>
  document.addEventListener("DOMContentLoaded", function() {
    var ctx = document.getElementById('myChart').getContext('2d');
    var myChart = new Chart(ctx, {
      type: 'line',
      data: {
        labels: [1, 2, 3, 4], // x values
        datasets: [{
          label: 'Simple Line Chart',
          data: [2, 4, 6, 8], // y values
          borderColor: 'rgba(75, 192, 192, 1)',
          backgroundColor: 'rgba(75, 192, 192, 0.2)',
          borderWidth: 1
        }]
      },
      options: {
        scales: {
          xAxes: [{
            scaleLabel: {
              display: true,
              labelString: 'X values'
            }
          }],
          yAxes: [{
            scaleLabel: {
              display: true,
              labelString: 'Y values'
            },
            ticks: {
              beginAtZero: true
            }
          }]
        }
      }
    });
  });
</script>

```