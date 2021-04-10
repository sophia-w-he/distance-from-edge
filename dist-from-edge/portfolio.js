"use strict";
(function() {

  window.addEventListener('load', init);

  function init() {

    let myimage = document.getElementById('myimage');
    console.log(myimage)

    myimage.onclick = function(e) {
      var ratioX = e.target.naturalWidth / e.target.offsetWidth;
      var ratioY = e.target.naturalHeight / e.target.offsetHeight;

      var domX = e.x + window.pageXOffset - e.target.offsetLeft;
      var domY = e.y + window.pageYOffset - e.target.offsetTop;

      var imgX = Math.round(domX * ratioX);
      var imgY = Math.round(domY * ratioY);

      console.log(imgX, imgY);
      postData(imgX, imgY);
    };



  }

  /*function changeImage(imgName, imgX, imgY) {
    jQuery.get('/dist-from-edge/distedge.py', function(data) {
    image = document.getElementById('imgDisp');
    image.src = data;
  }*/

  function postData(input1, input2) {
    /*$.ajax({
        type: "POST",
        url: "/dist-from-edge/distedge.py",
        data: { param: input },
        //dataType: "text",
        success: callbackFunc
    });*/
    $.ajax({
      type: 'GET',
      url: "http://127.0.0.1:5000/get_result/" + input1 + "/" + input2,
      data: {field: input1}, //passing some input here
      dataType: "text",
      success: function(response) {
          refreshImage('result-image', "dist-from-edge-test.jpg");
          alert("The minimum distance from your point to the closest edge is " + response);
          /*let resultImage = document.getElementById('result-image');
          console.log(resultImage)
          resultImage.src = "dist-from-edge-test.jpg";*/
      }
    }).done(function(data){
        console.log(data);
    }).fail(function() {
    alert( "Error" );
    })
    .always(function() {
      alert( "Finished! See results below." );
    });
  }

  function callbackFunc(response) {
      // do something with the response
      console.log(response);
  }

  function refreshImage(imgElement, imgURL){
    // create a new timestamp
    var timestamp = new Date().getTime();
    var el = document.getElementById(imgElement);
    var queryString = "?t=" + timestamp;
    el.src = imgURL + queryString;
  }

  //postData('data to process');

})();
