require.config({
        paths: {
          cytoscape: 'https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.2.8/cytoscape',
        }
      })

      require(['cytoscape'], function(cytoscape){
        $("#cyDiv").remove();
        element.append("<div id='cyDiv'></div>");
        $("#cyDiv").height("300px");

          var cy = cytoscape({
            container : document.getElementById('cyDiv'),
            elements  : %s,
            style     : %s,
            layout    : %s,
            maxZoom   : 10,
            minZoom   : 0.1,
          });
      });
