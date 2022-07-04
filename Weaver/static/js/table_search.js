$(document).ready(function(){
  $("#table_search-input").on("keyup paste", function(e) {
    var element = $(this);
    var event = e;
    setTimeout(function() {
        if (event.which == 27) { // Esc
            $(this).val("");
        }
        var value = element.val().toLowerCase();
        $(".table-search-target tbody tr").filter(function() {
            if($(this).find('.table-search-search_on').first().attr('data-search').toLowerCase().indexOf(value) > -1){
                $(this).removeClass('table_search-hide');
            } else {
                $(this).addClass('table_search-hide');
            }
        });
    }, 100);
  });
  $("#table_search-clear").click(function() {
    $(".table-search-target tbody tr").removeClass('table_search-hide');
    $("#table_search-input").val('');
  });
});