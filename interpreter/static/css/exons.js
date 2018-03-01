$(document).ready(function(){
    $("#1").click(function(){
        $("#1").fadeOut();
        $r='exon 1';
        $if(r.addEventListener){
            r.addEventListener('change',handleUserAnswer);
        }
        $else if(r.attachEvent){
            r.attachEvent('onchange',handleUserAnswer);
        }
      });

//     $("#2").click(function(){
//         $("#2").fadeOut();
//       });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
//     $("#3").click(function(){
//         $("#3").fadeOut();
//     });
// });
function handleUserAnswer(e){
    var el = e.target || e.srcElement;
    document.getElementById('userAnswer').value = document.getElementById('userAnswer').value.concat(el.value);
}

        if(r.addEventListener){
            r.addEventListener('change',handleUserAnswer);
        }else if(r.attachEvent){
            r.attachEvent('onchange',handleUserAnswer);
        }else{
            r.onChange = handleUserAnswer;
        }
    }
}