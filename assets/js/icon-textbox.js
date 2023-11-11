// This code waits until the DOM is fully loaded before running
document.addEventListener('DOMContentLoaded', (event) => {
  // Function to toggle the text box display
  function toggleTextBox(id) {
    var textBox = document.getElementById(id);
    if (textBox.style.display === "none") {
      textBox.style.display = "block";
    } else {
      textBox.style.display = "none";
    }
  }

  // Attach event listeners to your icons
  document.querySelectorAll('.icon').forEach(item => {
    item.addEventListener('click', event => {
      event.preventDefault();
      const textBoxId = event.target.getAttribute('href').substring(1);
      toggleTextBox(textBoxId);
    });
  });
});
